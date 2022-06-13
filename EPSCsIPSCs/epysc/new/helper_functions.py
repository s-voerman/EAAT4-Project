import numpy as np
from scipy.signal import (firwin, lfilter, butter, deconvolve, filtfilt,
                          savgol_filter, sosfilt, detrend, resample,
                          medfilt, wiener)
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from neo.rawio import AxonRawIO

# some helper functions and basic stuff for the processing of ABf files.


def detrend_data(data):
    return np.asarray(detrend(data))


def wiener_filter(data, order):
    return wiener(data, order)


def resample_data(data, resampling_rate, sr):
    new_total_size = int((len(data)/sr) * resampling_rate)
    data = resample(data, new_total_size)
    return data, resampling_rate


def exp_func(x, a, tau, c):
    return a * -np.exp(-x/tau) + c


def decay_fit_func(x, y, expected_values):
    return curve_fit(exp_func, x, y, expected_values)


def gaussian(data_decon, gaussian_size):
    gaussian_filtered = gaussian_filter(data_decon, gaussian_size)
    return gaussian_filtered


def filtering(data, cutoff, btype, sr, order=5):
    nyq = sr * 0.5
    cutoff /= nyq
    b, a = butter(order, cutoff, btype=btype, analog=False)
    y = filtfilt(b, a, data)  # goes to self.detection_data

    cutoff = 500 / nyq
    b, a = butter(order, cutoff, btype=btype, analog=False)
    orig_y = filtfilt(b, a, data)  # goes to self.filtered_data
    return y, orig_y


def read_abf(file_loc):
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header["signal_channels"][0][2])
    data = a.get_analogsignal_chunk()
    wrong_gain = a.header["signal_channels"][0][5]
    gain_correction = a._axon_info["fTelegraphAdditGain"][0]
    #correct_gain = wrong_gain * gain_correction
    #data = data * correct_gain
    return data[:, 0], sr

def read_abf(file_loc):
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header["signal_channels"][0][2])
    data = a.get_analogsignal_chunk() * a.header["signal_channels"][0][5]
    return data[:, 0], sr
