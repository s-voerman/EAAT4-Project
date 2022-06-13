import numpy as np
import os
import matplotlib.pyplot as plt
import datetime
from scipy.signal import (firwin, lfilter, butter, deconvolve, filtfilt,
                          savgol_filter, sosfilt, detrend, resample,
                          medfilt, wiener)
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from statsmodels.robust import mad
import pandas as pd
from neo.rawio import axonrawio
import helper_functions as hf
from process_event import analyze_event

plt.ion()

# implement init template creation
# for curve fit: fix x to [0,1,2...]
# also add p0 for curve_fit
# TODO: implement fitting, then .37 determination. should work better.
# note: Remove cell 030718 2-19.ABF, cause it died.


class process_ABF(object):
    def __init__(self, file_loc, output_loc, plot=False):
        # parameters
        self.file_loc = file_loc
        self.output_loc = output_loc
        self.plot = plot
        self.sr = 1
        self.original_sr = 1
        self.multiplier = 1 #was 3.
        self.threshold = None
        self.file_length = 0

        # data stores
        # to be clear: all detection is done on detection_data; all fitting etc
        # is done on filtered_data.
        self.original_data = None  # contains untouched data
        self.filtered_data = None  # contains filtered data
        self.detection_data = None  # contains detrended, resampled data
        self.deconvoluted_data = None  # contains results of deconvolution

    def create_template_vars(self):
        # old: r .0024, d .008
        # old 240718: r: 0.003; t = 0.01; a = -15
        # 011019: r: 0.0025, d: 0.009, -15
        self.rise = 0.0025
        self.decay = 0.009
        self.amplitude = -45

    def get_detection_vars(self):
        # event detection 
        # ideally 20 ms
        self.min_distance = 0.010 * self.sr

    def create_template(self, rise=None, decay=None, amplitude=None, plot=False):
        """ to properly deconvolve the data, we need a reliable
        template. note that the length of the template should be
        relatively long, as short templates induce oscillations """
        self.create_template_vars()
        if rise is None:
            rise = self.rise
        if decay is None:
            decay = self.decay
        if amplitude is None:
            amplitude = self.amplitude
        rise = rise * self.sr
        decay = decay * self.sr
        amp_normalization = (decay / rise) ** (rise / (rise - decay))
        template_size = round(0.05 * self.sr)
        g_x = np.arange(0, template_size)

        normalized_amplitude = amplitude / amp_normalization
        g_y = [normalized_amplitude * (-np.exp(-t/rise) + np.exp(-t/decay))
               if t > 0 else 0 for t in g_x]
        if plot:
            fig, ax = plt.subplots(figsize=(2, 2))
            ax.plot(g_x, g_y)
            plt.tight_layout()
        return g_y, g_x

    def deconvolution(self, data, template, wiener_flag=True,
                      wiener_order=255):
        if wiener_flag:
            data = hf.wiener(data, wiener_order)
        minute_in_samples = 20 * self.sr
        data_decon = np.array([])
        minute_windows = np.arange(0, len(data), minute_in_samples)
        for i, x in enumerate(minute_windows):
            print("window", i)
            start = x
            if i == len(minute_windows)-1:
                end = len(data)
            else:
                end = x + minute_in_samples
            dd = deconvolve(data[start:end], np.asarray(template[2:]))[0]
            data_decon = np.concatenate([data_decon, dd])
        dd = None
        data_decon[np.where(data_decon > 1.5)] = 0
        data_decon[np.where(data_decon < -1.5)] = 0
        return data_decon, data

    def get_threshold(self, data_decon):
        thresholds = []
        for x in np.arange(0, self.file_length, 3):
            x = x * self.sr
            step = 3 * self.sr
            m = np.median(data_decon[x:x+step])
            s = np.std(data_decon[x:x+step])
            threshold = m + self.multiplier*s
            thresholds.append(np.abs(threshold))
        threshold = np.min(thresholds)
        return threshold

    def propagate_threshold(self, data_decon):
        self.get_detection_vars()
        threshold = self.get_threshold(data_decon)
        # above_thresh = list(np.where(data_decon > threshold)[0])
        above_thresh = []
        for i, e in enumerate(data_decon[1:-1]):
            i = i + 1
            if (((e - np.mean(data_decon[i-50:i])) > threshold)
                and ((e - np.mean(data_decon[i+1:i+51])) > threshold)):
                above_thresh.append(i)

        indices = [i + 1 for (x, y, i) in zip(above_thresh, above_thresh[1:],
                                              range(len(above_thresh))) if
                   self.min_distance < abs(x - y)]
        events = [above_thresh[s:e] for s, e in zip([0] + indices, indices +
                                                    [len(above_thresh)])]

        output_events = np.array([])
        if len(events) == 1:
            if len(events[0]) == 0:
                return None, None
        for e in events:
            e = np.array(e).astype(int)
            m = np.max(data_decon[e])
            event_peak = e[np.where(data_decon[e] == m)[0]]
            output_events = np.append(output_events, event_peak[0])

        delta = len(self.detection_data) - len(data_decon)
        single_delta = delta/6
        print(output_events[:-3])
        output_events = np.array([x + (single_delta * int(x / (20*self.sr))) for x in output_events])
        print(output_events[:-3])
        # output_events = np.array([x if x < len(data_decon)/2 else x + delta / 2 for x in output_events])
        return events, output_events, threshold

    def quick_plot(self, threshold, peaks=None, plot_threshold=True):
        fig, ax = plt.subplots(3, 1, sharex=True)
        for i, y in enumerate([self.filtered_data,
                               self.detection_data,
                               self.deconvoluted_data]):
            ax[i].plot(np.linspace(0, self.file_length, len(y)), y, linewidth=.5)
        ax[0].set_xlim(0, 1)
        ax[0].set_title(self.file_loc)
        if plot_threshold:
            ax[2].axhline(threshold, color="red", linestyle=":")
        if peaks.any():
            peaks = peaks / self.sr
            peaks = peaks * self.original_sr
            peaks = peaks.astype(int)
            ax[0].scatter(np.linspace(0, self.file_length, len(self.filtered_data))[peaks],
                          self.filtered_data[peaks], color="green", s=5)
        plt.tight_layout()
        return fig, ax

    def output_plots(self, output_dict):
        pdo = pd.DataFrame()
        for x in output_dict:
            pdo = pdo.append(output_dict[x], ignore_index=True)

        filename = self.file_loc.split("/")[-1][:-4]
        print(len(self.filtered_data), self.original_sr)
        for x in np.arange(0, len(self.filtered_data), self.original_sr):
            start = int(x / self.original_sr)
            end = start + 1
            z = pdo.loc[pdo.time.between(start, end, inclusive=False)]
            fig, ax = plt.subplots(figsize=(15,5))
            ax.plot(np.linspace(start, end, self.original_sr),
                    self.filtered_data[x:x+self.original_sr], linewidth=.4)
            ax.set_xlim(start, end)
            ax.set_ylim(ax.get_ylim()[0]-10, ax.get_ylim()[1]+10)

            # select reasonable events
            # z = z.loc[np.abs(z.amplitude) > z.rms]
            # z = z.loc[z.tau > 0.0001]
            # z = z.loc[z.rise > 0.00004]
            # z = z.loc[z.fit_tau < 10000]  # 10000
            # z = z.loc[z.fit_tau > 0]
            # z = z.loc[z.tau > z.rise * 0.9]
            # z = z.loc[np.abs(z.tau_pams) > 100]

            if z.size > 0:
                for i, event in z.iterrows():
                    color = "orange"
                    if np.abs(event.amplitude) > 2 * event.rms:
                        color = "green"
                    xs = event["time"]
                    xp = event["rise"] + xs
                    xt = event["tau"] + xp
                    yb = event["baseline"]
                    yp = event["amplitude"]
                    yt = yb + (np.exp(-1)*yp)
                    fit_x = np.array((event["fit_x"]) / self.original_sr) + xp
                    fit_y = event["fit_y"]
                    ax.scatter(xs, yb, color=color, s=20)
                    ax.scatter(xp, yb + yp, color=color, s=20)
                    ax.scatter(xt, yt, color=color, s=20)
                    ax.plot(fit_x, fit_y, linewidth=1, color="red")
            ax.set_title(filename)
            base = "./pics/" + self.output_loc + "/"
            if not os.path.exists(base+filename):
                os.makedirs(base+filename)
            fig.savefig(base + filename + "/" + filename + "_" +
                        str(int(x/self.original_sr)).zfill(3) + ".png")
            plt.close()
        return


    def run_all(self):
        self.original_data, self.sr = hf.read_abf(self.file_loc)
        self.original_sr = int(self.sr)

        print(1)
        self.detection_data, self.filtered_data = hf.filtering(self.original_data, 200, "lowpass", self.sr)
        print(2)
        self.file_length = int(len(self.filtered_data) / self.sr)
        print(3)
        self.detection_data, self.sr = hf.resample_data(self.detection_data,
                                                        2500, self.original_sr)
        print(4)

        self.detection_data = hf.detrend_data(self.detection_data)
        print(5)
        template, template_y = self.create_template()
        print(6)
        self.deconvoluted_data, self.detection_data = self.deconvolution(self.detection_data,
                                                                         template,
                                                                         wiener_flag=True,
                                                                         wiener_order=13)
        print(7)
        self.deconvoluted_data = np.abs(self.deconvoluted_data)
        print(8)
        self.deconvoluted_data = hf.gaussian(self.deconvoluted_data, 3)
        print(9)
        events, outevs, threshold = self.propagate_threshold(self.deconvoluted_data)

        if not events:
            return None
        if self.plot:
            fig, ax = self.quick_plot(threshold, outevs)
        else:
            fig, ax = None, None

        output_dict_of_dicts = {}
        print("analysing events")
        for i, x in enumerate(outevs):
            analyze = analyze_event(self.plot, fig, ax, x, self.filtered_data, self.original_sr, self.sr)
            analyze_output = analyze.process_event()
            output_dict_of_dicts[i] = analyze_output

        #self.output_plots(output_dict_of_dicts)
        #temp
        if self.plot:
            return output_dict_of_dicts, self.filtered_data, self.deconvoluted_data, self.detection_data, fig, ax
        else:
            return output_dict_of_dicts


