import os
from scipy.stats import spearmanr
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from neo.rawio import axonrawio
import helper_functions as hf
from scipy.signal import savgol_filter

'''
Based on:
    Counting channels: a tutorial guide on channel fluctuation analysis by
    Alvarez, Gonzales and Latorre (2002)
    Studying properties of neurotransmitter receptors by non-stationary noise
    analysis of spontaneous postsynaptic currents and agonist-evoked responses
    in outside-out patches by Hartveit and Veruki, 2007

Implement changes upstream:
    Make rise time 10-90%. Seems a bit more reliable possibly? Maybe also add a
    quick fit to the rise time? Seems like it could be useful, but a bit wonky
    maybe, when you get double events or elongated event peaks (you know what i
    mean).

Basic setup:
    Find locations of all events, based on xlsx data from EPSC analysis.
    For each file, plot the per-event amplitude, decay, rise versus time as x.
    This allows you to detect (hopefully) patterns in the change of event
    characteristics, which tells you something about the events. Save these
    with appropriate title and filename in a folder to look at. 
    Extract all events, with appropriate tau length.
    Remove events where the events are too close to others (i.e. unstable pre-
    and post-event).
    Pad all events at the start and end to avoid problems with event lengths.
    Calculate the mean event.
    Normalize event to mean event, calculate variance per event and per
    timepoint. Preferably use bins, e.g. .5 ms long. 
    Do other stuff, presumable?

Test-data:
    ./data/190718/mepsc/190718 5-43.ABF for wt
    ./data/240718/mepsc/240718 6-37.ABF for mutant
'''


class NSFA(object):
    def __init__(self, file_loc, event_xlsx, output_folder):
        self.file_loc = file_loc
        self.sr = None
        self.event_data = None
        self.event_xlsx = event_xlsx
        self.abf_data = None
        self.extracted_events = {}
        self.event_metadata = None
        self.output_folder = output_folder

    def get_extraction_parameters(self):
        self.padding = 0.002 * self.sr

    def generate_summary_graph(self, data):
        data.tau = data.tau * 1000
        data.rise = data.rise * 1000
        fig, ax = plt.subplots(2, 3, figsize=(10, 6))
        p_values = []
        for i, y in enumerate(["amplitude", "rise", "tau"]):
            sns.regplot(x="time", y=y, data=data, scatter_kws={"s": 3},
                        ax=ax[0, i])
            rho, p = spearmanr(data["time"], data[y])
            p_values.append(p)
            ax[0, i].text(0.1, 0.1, "rho = %.3f" % rho, transform=ax[0, i].transAxes)
            ax[0, i].text(0.1, 0.2, "p = %.3f" % p, transform=ax[0, i].transAxes)
        for i, c in enumerate(zip(["amplitude", "rise", "tau"], ["rise", "tau",
                                                                 "amplitude"])):
            sns.regplot(x=c[0], y=c[1], data=data, scatter_kws={"s": 3},
                        ax=ax[1, i])
            rho, p = spearmanr(data[c[0]], data[c[1]])
            p_values.append(p)
            ax[1, i].text(0.1, 0.1, "rho = %.3f" % rho,
                          transform=ax[1,i].transAxes)
            ax[1, i].text(0.1, 0.2, "p = %.3f" % p,
                          transform=ax[1,i].transAxes)

        fig.suptitle(x)
        plt.tight_layout()
        plt.subplots_adjust(top=0.92)
        filename = self.file_loc.split("/")[-1][:-4]
        if not os.path.exists(self.output_folder + "summary_graphs/"):
            os.makedirs(self.output_folder + "summary_graphs/")
        fig.savefig(self.output_folder + "summary_graphs/" + filename + ".png")
        plt.close()
        return np.array(p_values)

    def get_event_data(self):
        data = pd.read_excel(self.event_xlsx)
        self.event_metadata = data.loc[data.id == self.file_loc]
        p_values = self.generate_summary_graph(self.event_metadata)
        if (p_values < 0.05).any():
            print("Warning: significant pattern detected")
        return

    def read_abf(self, filtering=False, resample=False):
        self.abf_data, self.sr = hf.read_abf(self.file_loc)
        if filtering:
            self.abf_data, _ = hf.filtering(self.abf_data, 750, "low", self.sr)
        if resample:
            self.abf_data, self.sr = hf.resample_data(self.abf_data, 2500,
                                                      self.sr)
        return

    def cut_out_events(self):
        self.get_extraction_parameters()
        for i, r in self.event_metadata.iterrows():
            event_info = {}
            start_time = r["time"]
            total_length = r["rise"] + 2 * r["tau"]
            end_time = start_time + total_length
            print(start_time, end_time)
            start_time_s = int(start_time * self.sr)
            end_time_s = int(end_time * self.sr)
            print(start_time_s, end_time_s)
            if self.padding:
                start_time_s -= int(self.padding)
            baseline = np.mean(self.abf_data[start_time_s:start_time_s+int(self.padding)])

            self.extracted_events[i] = (self.abf_data[start_time_s:end_time_s]
                                        - baseline)
        return

    def pad_events(self):
        max_len = 0
        for x in self.extracted_events:
            if len(self.extracted_events[x]) > max_len:
                max_len = len(self.extracted_events[x])
        for x in self.extracted_events:
            self.extracted_events[x] = (np.pad(self.extracted_events[x].astype(float),
                                               [0, max_len-len(self.extracted_events[x])],
                                               "constant", constant_values=np.nan))
        return

    def line_up_graph(self):
        fig, ax = plt.subplots()
        for x in self.extracted_events:
            d = self.extracted_events[x]
            if self.sr == 50000:
                win = 501
            else:
                win = 25
            d = savgol_filter(d, win, 3, mode="mirror")
            min_loc = np.where(np.diff(d) == np.min(np.diff(d)))[0][0]
            x_vals = np.arange(-min_loc, (len(d)-min_loc))
            ax.plot(x_vals, d, linewidth=.4)
        return fig, ax

    def run(self):
        self.read_abf(resample=True)
        self.get_event_data()
        # self.get_extraction_parameters()
        # self.cut_out_events()
        # self.pad_events()
        fig, ax = self.line_up_graph()
        return self.abf_data, self.extracted_events, fig, ax

if __name__ == "__main__":
    file_loc = "./data/190718/mepsc/190718 3-23.ABF"
    metadata = "./output/final_270718/cumulative_output.xlsx"
    output_folder = "./NSFA/050818/"
    d = pd.read_excel(metadata)
    file_locs = list(d.id.unique())
    for x in file_locs:
        file_loc = x
        analysis = NSFA(file_loc, metadata, output_folder)
        data, events, fig, ax = analysis.run()
