import numpy as np
import os
import matplotlib.pyplot as plt
import helper_functions as hf
import pandas as pd

class analyze_event(object):
    """
    Receive a single event (i.e. the peak determined via threshold search
    after deconvolution. Find the peak, baseline, and rise. Fit the decay with
    high temporal resolution, and find the tau. Fancy boyes.

    input vars:
        filtered_data = array containing the filtered but otherwise unprocessed
        data.
        event = timepoint of event in samples
    output:
        amplitude, baseline, rise, tau, start_time = all floats
        fit parameters = array
    """

    def __init__(self, plot, fig, ax, event, filtered_data, data_sr, detection_sr):
        self.plot = plot
        self.fig = fig
        self.ax = ax
        self.detection_sr = detection_sr
        self.data = filtered_data
        self.data_sr = data_sr
        self.rms_segment_size = 2 * self.data_sr
        self.event = int((event / self.detection_sr) * self.data_sr)

    def generate_detection_parameters(self):
        # detection parameters
        self.peak_search = int(round(0.005 * self.data_sr))
        self.baseline_length = int(round(0.002 * self.data_sr))
        self.decay_threshold = 0.2
        self.decay_search = int(0.050 * self.data_sr)
        self.decay_extra_length = int(0.010 * self.data_sr)

    def find_peak(self):
        x = self.event
        i = np.arange(x-self.peak_search, x+self.peak_search).astype(int)
        peak_current = np.min(self.data[i])
        location_of_peak = i[np.where(self.data[i] == peak_current)].item()
        return peak_current, location_of_peak

    def find_baseline(self, location_of_peak):
        # find the start of the event, average current in the x ms before the
        # event, and return these values as well as the x_values for the
        # baseline_current.
        x = location_of_peak - 10
        while (self.data[x] - np.mean(self.data[x-200:x]) < 0):
            x = x - 1
        event_start = x

        start_current = self.data[x]
        baseline_range = np.arange(x-self.baseline_length, x)
        baseline_current = np.mean(self.data[baseline_range.astype(int)])
        return baseline_current, start_current, event_start, baseline_range

    def find_amplitude(self, peak_current, baseline_current, start_current):
        amplitude = peak_current - baseline_current
        exact_amplitude = peak_current - start_current
        return amplitude, exact_amplitude

    def fit_decay(self, location_of_peak, exact_amplitude, start_current,
                  baseline_current):
        decay_search_current = start_current + (self.decay_threshold *
                                                exact_amplitude)
        decay_search = (0.025 + (np.abs(exact_amplitude) / 1000)) * self.data_sr
        decay_search_x_values = np.arange(location_of_peak, location_of_peak+int(decay_search))

        try:
            location_of_end_decay = np.where(self.data[decay_search_x_values] >
                                             decay_search_current)[0][0]
        except IndexError:
            return None
        decay_end = decay_search_x_values[location_of_end_decay] + self.decay_extra_length
        decay_y_values = self.data[int(location_of_peak):int(decay_end)]
        decay_x_values = np.arange(len(decay_y_values))
        expected_values = [1,len(decay_x_values)/2, baseline_current]

        # fit with decay function
        try:
            fit, coeff = hf.decay_fit_func(decay_x_values, decay_y_values,
                                    expected_values)
        except (TypeError, RuntimeError) as e:
            return None

        high_def_x_values = np.arange(0, len(decay_y_values), 0.01)
        high_def_y_values = hf.exp_func(high_def_x_values, fit[0], fit[1],
                                        fit[2])
        decay_tau_threshold = start_current + (np.exp(-1) * exact_amplitude)

        # CALCULATE sum of squared deviations
        difference = decay_y_values - high_def_y_values[::100]
        difference = difference ** 2
        difference = np.sum(np.sqrt(difference)) / len(decay_y_values)

        try:
            tau = high_def_x_values[np.where(high_def_y_values > decay_tau_threshold)][0]
        except IndexError:
            return None
        return tau, fit, decay_tau_threshold, high_def_x_values, high_def_y_values, coeff, difference

    def determine_noise(self, location_of_peak):
        segment = np.arange(location_of_peak - self.rms_segment_size * 0.5,
                            location_of_peak + self.rms_segment_size * 0.5).astype(int)
        if segment[0] < 0:
            segment = segment[np.where(segment >= 0)]
        if segment[-1] > len(self.data):
            segment = segment[np.where(segment < len(self.data))]

        segment_data = self.data[segment]
        segment_data = hf.detrend(segment_data)
        rms_list = []
        steps = int(self.data_sr * 0.1)
        for x in np.arange(0, len(segment_data), steps):
            rms = np.sqrt(np.mean(segment_data[x:x+steps] ** 2))
            rms = rms * (2 * np.sqrt(2))
            rms_list.append(rms)

        # perc_50 = np.percentile(rms_list, 25)
        # rms_list = np.array(rms_list)
        # rms_list = rms_list[np.where(rms_list < perc_50)]
        return np.min(rms_list)

    def plot_event_parameters(self, output_dict):
        if output_dict["amplitude"] > output_dict["rms"] * -1.5:
            return
        elif output_dict["tau"] == None:
            return
        elif output_dict["tau"] < 0.001:
            return

        xs = output_dict["time"]
        xp = output_dict["rise"] + xs
        xt = output_dict["tau"] + xp
        yb = output_dict["baseline"]
        yp = output_dict["exact_amplitude"]
        fit_x = output_dict["fit_x"] / 50000
        fit_y = output_dict["fit_y"]
        yt = yb + (np.exp(-1)*yp)

        #ratio2_string = "pa/ms ratio = %.3f" % (output_dict["rise_pams"] / output_dict["tau_pams"])
        #tau_string = "tau pa/ms = %.3f" % round(output_dict["tau_pams"], 3)
        #rise_string = "rise pa/ms = %.3f" % round(output_dict["rise_pams"], 3)
        #ratio_string = "ratio = %.3f" % round(output_dict["rise_tau_ratio"], 3)
        #self.ax[0].text(xs, (yb+yp) + 23, rise_string, fontsize=7)
        #self.ax[0].text(xt, (yb+yp) + 26, tau_string, fontsize=7)
        #self.ax[0].text(xp, (yb+yp) + 30, ratio_string, fontsize=7)
        #self.ax[0].text(xp, (yb+yp) + 40, ratio2_string, fontsize=7)
        self.ax[0].scatter(xs, yb, color="orange")
        self.ax[0].scatter(xp, yb + yp, color="orange")
        self.ax[0].scatter(xt, yt, color="orange")
        self.ax[0].plot(np.array(fit_x) + xp, fit_y, linewidth=1, color="purple")

        return

    def generate_output_dict(self, amplitude, exact_amplitude,
                             baseline_current, location_of_peak,
                             event_start, start_current, tau, fit,
                             fit_x, fit_y, coeff, difference):
        output_dict = {}
        output_dict["time"] = event_start / self.data_sr
        output_dict["amplitude"] = amplitude
        output_dict["exact_amplitude"] = exact_amplitude
        output_dict["baseline"] = baseline_current
        output_dict["rise"] = (location_of_peak - event_start) / self.data_sr
        output_dict["rise_pams"] = (amplitude / ((location_of_peak - event_start) /
                                    self.data_sr))
        output_dict["fit_x"] = fit_x
        output_dict["fit_y"] = fit_y
        output_dict["fit_quality"] = difference
        for i, x in enumerate(["fit_coefficient", "fit_tau", "fit_constant"]):
            try:
                output_dict[x] = fit[i]
            except TypeError:
                output_dict[x] = None

        for i, x in enumerate(["fit_coefficient_variance",
                               "fit_tau_variance",
                               "fit_constant_variance"]):
            try:
                output_dict[x] = coeff[i,i]
            except TypeError:
                output_dict[x] = None

        if tau:
            output_dict["tau_orig"] = tau
            output_dict["tau"] = tau / self.data_sr
            output_dict["tau_pams"] = amplitude / (tau / self.data_sr)
            output_dict["rise_tau_ratio"] = output_dict["rise"] / output_dict["tau"]
        else:
            output_dict["tau"] = None
            output_dict["tau_orig"] = None
            output_dict["tau_pams"] = None
            output_dict["rise_tau_ratio"] = None
        output_dict["rms"] = self.determine_noise(location_of_peak)
        return output_dict

    def process_event(self):
        self.generate_detection_parameters()
        peak_current, location_of_peak = self.find_peak()
        baseline_current, start_current, event_start, baseline_range = self.find_baseline(location_of_peak)

        amplitude, exact_amplitude = self.find_amplitude(peak_current,
                                                         baseline_current,
                                                         start_current)
        o = self.fit_decay(location_of_peak, exact_amplitude, start_current,
                           baseline_current)
        try:
            tau, fit, decay_tau_threshold, fit_x, fit_y, coeff, difference = o
        except TypeError:
            tau, fit, decay_tau_threshold, fit_x, fit_y, coeff, difference = [None] * 7

        output_dict = self.generate_output_dict(amplitude, exact_amplitude,
                                                baseline_current,
                                                location_of_peak,
                                                event_start, start_current,
                                                tau, fit, fit_x, fit_y,
                                                coeff, difference)
        if self.plot:
            self.plot_event_parameters(output_dict)
        return output_dict
