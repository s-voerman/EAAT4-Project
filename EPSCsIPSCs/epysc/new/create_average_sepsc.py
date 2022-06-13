import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import helper_functions as hf
from scipy.stats import sem
from scipy.signal import (firwin, lfilter, butter, deconvolve, filtfilt,
                          savgol_filter, sosfilt, detrend, resample,
                          medfilt, wiener)
from neo.io import AxonIO

"""
Find all events, line em up and extract em. Then create an average trace.
"""

data = pd.read_excel("./all_compensated_data_280618.xlsx", sheetname=1)
prefix = "./data/sEPSC/"
data.loc[:, "id"] = [prefix + str(x).zfill(6) + " " + y + ".ABF" for
                     x, y in zip(data.date, data.cell)]
data = data[["id", "group", "animal", "Time (ms)", "Amplitude", "Rise (ms)",
             "Decay (ms)", "Baseline", "HalfWidth", "date"]]
c = list(data.columns)
c[3] = "time"
c[5] = "rise"
c[6] = "decay"
data.columns = c
data.loc[:, "time"] /= 1000
data.loc[:, "decay"] /= 1000
data.loc[:, "rise"] /= 1000


data.loc[:, "sr"] = [20000 if str(x) != "240117" else 10000 for x in data.date]


def detrend_data(data):
    return np.asarray(detrend(data))

def get_abf(source):
    block = AxonIO(source).read_block()
    data = np.array(block.segments[0].analogsignals[0])
    return data

def filtering(data, cutoff, btype, sr, order=5):
    nyq = sr * 0.5
    cutoff /= nyq
    b, a = butter(order, cutoff, btype=btype, analog=False)
    y = filtfilt(b, a, data)  # goes to self.detection_data

    return y

def resample_data(data, resampling_rate, sr):
    new_total_size = int((len(data)/sr) * resampling_rate)
    data = resample(data, new_total_size)
    return data, resampling_rate

def get_original_trace(data):
    rsr = 5000.0
    full_length = 0.02 * rsr
    padding = 0.01 * rsr
    out = []
    # data = data.loc[(data.amplitude < np.percentile(data.amplitude, 95)) &
    #                 (data.amplitude > np.percentile(data.amplitude, 5)) &
    #                 (data.tau < np.percentile(data.tau, 95)) &
    #                 (data.tau > np.percentile(data.tau, 5))]
    for d in data.id.unique():
        sr = data.loc[data.id == d, "sr"][0]
        trace = get_abf(d).astype(float)
        trace = np.hstack(np.array(trace))
        trace, _ = resample_data(trace, 5000, sr)
        trace = filtering(trace, 500, "lowpass", 5000)
        # trace, _ = filtering(trace, 5, "highpass", 5000)
        trace = detrend_data(trace)
        cell_data = data.loc[data.id == d]
        cell_data = cell_data.loc[cell_data.time < 119.80]
        cell_data.loc[:, "time"] *= rsr
        cell_data.loc[:, "rise"] *= rsr
        for i, r in cell_data.iterrows():
            peak = r.time - r.rise
            start = int(peak - padding)
            end = int(peak + full_length)
            t = np.array(trace[start:end])[0:150]
            if len(t) != 150:
                print(r, len(t))
            t = t/r.Amplitude
            t = np.array(t).astype(float)
            out.append(t)
    out = np.array(out)
    return out

def cleanup(d):
    clean = []
    for i in d:
        if len(i) != 150:
            continue
        i = np.array(i)
        avg = i[40:50].argmax()
        avg = 30 + avg
        i = i[avg - 20:]
        i = np.pad(np.array(i), (0, 150-len(i)), "constant", constant_values=(np.nan,
            np.nan))
        out = np.array(i)
        avg = np.nanmean(out)
        out = out - avg
        # out = out / peak
        # out = out * -1
        if np.isinf(out).any():
            continue
        clean.append(out[0:150])
    clean = np.array(clean)
    return clean

def cleanup(d):
    clean = []
    for i in d:
        if len(i) != 150:
            continue
        i = np.array(i)
        avg = i[40:50].argmax()
        avg = 30 + avg
        i = i[avg - 20:]
        i = np.pad(np.array(i), (0, 150-len(i)), "constant", constant_values=(np.nan,
            np.nan))
        out = np.array(i)
        avg = np.nanmean(out)
        out = out - avg
        # out = out / peak
        # out = out * -1
        if np.isinf(out).any():
            continue
        clean.append(out[0:150])
    clean = np.array(clean)
    return clean

def make_fig(data, dataset, fig=None, ax=None):
    if dataset == "wt":
        c = "black"
    else:
        c = "#9d1a20"
    if not fig:
        fig, ax = plt.subplots()
    f = np.nanmean(data, axis=0)
    minimum = np.min(f[~np.isnan(f)])
    f = (f / minimum)*-1
    f_sem = sem(data, axis=0)
    f_above = f + f_sem
    f_below = f - f_sem
    ax.plot(np.linspace(-10, 20, len(f)), f, color=c)
    ax.fill_between(np.linspace(-10, 20, len(f)), f_above, f_below,
                    color=c, alpha=0.5)
    return fig, ax

def find_prototype(data):
    out = {}
    for g in data.geno.unique():
        d = data.loc[data.geno == g]
        mean_amp = d.amplitude.mean()
        std_amp = d.amplitude.std()
        mean_tau = d.tau.mean()
        std_tau = d.tau.std()
        mean_rise = d.rise.mean()
        std_rise = d.rise.std()

        avg_dict = {}
        for v in zip(["amp", "tau", "rise"],
                     [mean_amp, mean_tau, mean_rise],
                     [std_amp, std_tau, std_rise]):
            bottom = v[1] - 0.218 * v[2]
            top = v[1] + 0.218 * v[2]
            avg_dict[v[0]] = [bottom, top]

        s = d.loc[(d.amplitude > avg_dict["amp"][0])
                  & (d.amplitude < avg_dict["amp"][1])]

        s = s.loc[(s.tau > avg_dict["tau"][0])
                  & (s.tau < avg_dict["tau"][1])]

        s = s.loc[(s.rise > avg_dict["rise"][0])
                  & (s.rise < avg_dict["rise"][1])]
        out[g] = s
    return out



