import matplotlib as mpl
import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import helper_functions as hf
import pandas as pd
from scipy.stats import levene, ttest_ind
from scipy.stats import linregress
from scipy.stats import tstd
from scipy.stats import spearmanr

plt.ion()
pd.options.display.max_columns = 20
pd.options.display.width = 160
pd.options.display.max_rows = 130

class process_single_cell(object):
    def __init__(self, file_loc, output_loc, cell, group, pre, post, row):
        self.file_loc = file_loc
        self.meta_row = row
        self.output_loc = output_loc
        self.cell = cell
        self.group = group
        self.pre = pre
        self.post = post
        # timepoints
        self.epsc1 = 0.1
        self.epsc2 = 0.2
        self.seal = [0.7, 0.8]
        self.sweep_size = 1.0
        return

    def process_file(self, data, sr, prepost):
        output = {}
        sweep_size = int(sr * self.sweep_size)
        windows = np.arange(0, len(data), sweep_size)
        output_pd = pd.DataFrame()
        print(len(data), sr, len(data)/sr)
        for i, x in enumerate(windows):
            output = {}
            d = data[x:x+sweep_size]
            output["pre_post"] = prepost
            if prepost == "post":
                output["time"] = i + 15
            else:
                output["time"] = i
            for epsc, start in zip(["EPSC1", "EPSC2"], [int(self.epsc1*sr),
                                                        int(self.epsc2*sr)]):
                baseline = np.mean(d[start-400:start-10])
                peak = np.min(d[start+1:start+400])
                amplitude = peak - baseline
                for sv, v in zip(["baseline", "peak", "amplitude"],
                                 [baseline, peak, amplitude]):
                    output["_".join([epsc, sv])] = v

            seal_baseline = np.mean(d[int(0.6*sr):int(0.695*sr)])
            seal_cap = np.min(d[int(0.69*sr):int(0.71*sr)]) - seal_baseline
            seal_input = np.mean(d[int(0.77*sr):int(0.795*sr)]) - seal_baseline
            output["r_cap"] = (-10e-3 / (seal_cap*1e-12)) / 1e6
            output["r_input"] = (-10e-3 / (seal_input*1e-12)) / 1e6
            output["r_membrane"] = output["r_input"] - output["r_cap"]
            output["holding"] = baseline
            output_pd = output_pd.append(output, ignore_index=True)
        output_pd.loc[:, "cell"] = self.cell
        output_pd.loc[:, "group"] = self.group
        return output_pd

    def run(self):
        pre_data = np.array([])
        for f in self.pre:
            filename = self.cell + "-" + f + ".ABF"
            data, sr = hf.read_abf(os.path.join(self.file_loc, filename))

            pre_data = np.concatenate([pre_data, data])

        post_data = np.array([])
        for f in self.post:
            filename = self.cell + "-" + f + ".ABF"
            data, sr = hf.read_abf(os.path.join(self.file_loc, filename))
            post_data = np.concatenate([post_data, data])

        sr = 50000
        pre_data = self.process_file(pre_data, sr, "pre")
        pre_data.loc[:, "time"] = [x - (np.max(pre_data.time)) for x in pre_data.time]
        pre_data.loc[:, "time"] = [int(x/3)-1 for x in pre_data.time]
        post_data = self.process_file(post_data, sr, "post")
        post_data.loc[:, "time"] = [int(x/3) for x in post_data.time]
        output_data = pd.concat([pre_data, post_data])
        #normalize
        output_data.loc[:, "PPR"] = output_data.EPSC2_amplitude / output_data.EPSC1_amplitude
        for i, v in enumerate(["EPSC1_amplitude", "EPSC2_amplitude",
                               "r_cap", "r_membrane", "PPR", "holding"]):
            output_data.loc[:, "norm_"+v] = (output_data.loc[:, v] /
                                            np.mean(output_data.loc[(output_data.time > -6) & (output_data.time < 0), v]))
        return output_data
    
def overview_graph(data, output_loc):
    #data = pd.DataFrame(data.groupby(["cell", "group", "time"]).mean().reset_index())
    data = pd.DataFrame(data.groupby(["cell", "time"]).mean().reset_index())
    for cell in data.cell.unique():
        d = data.loc[data.cell == cell]
        fig, ax = plt.subplots(6,2, sharex=True, figsize=(6,12))
        fig.suptitle(cell)
        labels = {"EPSC1_amplitude" : "EPSC1",
                  "EPSC2_amplitude" : "EPSC2",
                  "r_cap" : "Rseries",
                  "r_membrane" : "Rmembrane",
                  "PPR": "PPR",
                  "holding":"holding"}
        for i, v in enumerate(["EPSC1_amplitude", "EPSC2_amplitude",
                               "r_cap", "r_membrane", "PPR", "holding"]):
            ax[i,0].scatter(d.time, d[v], s=4)
            ax[i,0].set_xlabel("Time")
            ax[i,0].set_ylabel(labels[v])

            ax[i,1].scatter(d.time, d["norm_"+v], s=4)
            ax[i,1].axhline(0.8, lw=0.3)
            ax[i,1].axhline(1.2, lw=0.3)
            ax[i,1].axvline(25, lw=0.3, color="black")
            ax[i,1].axvline(-5, lw=0.3, color="black")
            ax[i,1].axvline(0, lw=0.3, color="black")
            ax[i,1].axvline(30, lw=0.3, color="black")
            ax[i,1].set_ylim(0,2)
            ax[i,1].set_xlabel("Time")
            ax[i,1].set_ylabel("norm" + labels[v])

        plt.tight_layout()
        plt.subplots_adjust(top=.92)
        fig.savefig(os.path.join(output_loc, "summary_"+cell+".png"))
        plt.close()
    return

def summary_graph(data, output_loc):
    #data = pd.DataFrame(data.groupby(["cell", "group", "time"]).mean().reset_index())
    data = pd.DataFrame(data.groupby(["cell", "time"]).mean().reset_index())
    for cell in data.cell.unique():
        d = data.loc[data.cell == cell]
        fig, ax = plt.subplots(6,2, sharex=True, figsize=(6,12))
        fig.suptitle(cell)
        labels = {"EPSC1_amplitude" : "EPSC1",
                  "EPSC2_amplitude" : "EPSC2",
                  "r_cap" : "Rseries",
                  "r_membrane" : "Rmembrane",
                  "PPR": "PPR",
                  "holding":"holding"}
        for i, v in enumerate(["EPSC1_amplitude", "EPSC2_amplitude",
                               "r_cap", "r_membrane", "PPR", "holding"]):
            ax[i,0].scatter(d.time, d[v], s=4)
            ax[i,0].set_xlabel("Time")
            ax[i,0].set_ylabel(labels[v])

            ax[i,1].scatter(d.time, d["norm_"+v], s=4)
            ax[i,1].axhline(0.8, lw=0.3)
            ax[i,1].axhline(1.2, lw=0.3)
            ax[i,1].axvline(25, lw=0.3, color="black")
            ax[i,1].axvline(-5, lw=0.3, color="black")
            ax[i,1].axvline(0, lw=0.3, color="black")
            ax[i,1].axvline(30, lw=0.3, color="black")
            ax[i,1].set_ylim(0,2)
            ax[i,1].set_xlabel("Time")
            ax[i,1].set_ylabel("norm" + labels[v])

        plt.tight_layout()
        plt.subplots_adjust(top=.92)
        fig.savefig(os.path.join(output_loc, "summary_"+cell+".png"))
        plt.close()
    return

def find_meta(meta_loc, data_loc):
    all_data = pd.DataFrame()
    meta = pd.read_excel(meta_loc)
    for i, row in meta.iterrows():
        group = row.group
        cell = " ".join([row["name"], str(row["cell"])])
        file_prefix = row.cell
        if isinstance(row.pre, str):
            pre = row.pre.split(", ")
        else:
            pre = [str(row.pre)]
        if isinstance(row.post, str):
            post = row.post.split(", ")
        else:
            post = [str(row.post)]
        p = process_single_cell(data_loc, "./output", cell, group, pre, post, row)
        output = p.run()

        all_data = pd.concat([all_data, output])
    return all_data

def convert_to_sasa_format(raw_data, good_cells=None):
    """
    for each file, we need to create a separate file. first sheet gives
    information on what animal etc it is. second sheet gives raw data (EPSC1,
    EPSC2, PPR, Rs, Ri, holding). Third gives the same, but per minute. Fourth
    gives normalized, fifth normalized per minute. Sounds good? Sounds good.

    fyi you get non-minutized data, and you can convert it to minutes by just
    groupbying with time.
    """
    out_norm = pd.DataFrame()
    out_raw = pd.DataFrame()
    out_PPR = pd.DataFrame()
    if good_cells:
        selector = good_cells
    else:
        selector = raw_data.cell.unique()
    for cell in selector:
        d = raw_data.loc[raw_data.cell == cell]
        info_file = {
            "cell": cell,
            "group": d.group.iloc[0],
        }
        info_cols = list(info_file.keys())
        info_sheet = d[info_cols]
        info_cols += ["time"]
        raw_full = [
            "EPSC1_amplitude",
            "EPSC2_amplitude",
            "holding",
            "PPR",
            "r_membrane",
            "r_cap"
        ]
        norm_full = ["norm_"+x for x in raw_full]
        raw_full = info_cols + raw_full
        norm_full = info_cols + norm_full
        raw_full_data = d[raw_full]
        raw_full_data_min = d.groupby(info_cols).mean().reset_index()[raw_full]
        norm_full_data = d[norm_full]
        norm_full_data_min = d.groupby(info_cols).mean().reset_index()[norm_full]

        # build excel file
        writer = pd.ExcelWriter("./output/100hz/"+cell+".xlsx")
        info_sheet.iloc[0].to_excel(writer, sheet_name="0_info")
        raw_full_data.to_excel(writer, sheet_name="1_raw")
        raw_full_data_min.to_excel(writer, sheet_name="2_raw_per_min")
        norm_full_data.to_excel(writer, sheet_name="3_norm")
        norm_full_data_min.to_excel(writer, sheet_name="4_norm_per_min")
        writer.save()

    writer = pd.ExcelWriter("./output/100hz/norm_epsc.xlsx")
    data = pd.DataFrame(raw_data)
    norm_epsc = pd.pivot_table(data, values="norm_EPSC1_amplitude",
                               index=["group", "cell"], columns="time").T
    norm_epsc.to_excel(writer)
    writer.save()

    writer = pd.ExcelWriter("./output/100hz/raw_epsc.xlsx")
    data = pd.DataFrame(raw_data)
    norm_epsc = pd.pivot_table(data, values="EPSC1_amplitude",
                               index=["group", "cell"], columns="time").T
    norm_epsc.to_excel(writer)
    writer.save()

    writer = pd.ExcelWriter("./output/100hz/raw_ppr.xlsx")
    data = pd.DataFrame(raw_data)
    norm_epsc = pd.pivot_table(data, values="PPR",
                               index=["group", "cell"], columns="time").T
    norm_epsc.to_excel(writer)
    writer.save()

    writer = pd.ExcelWriter("./output/100hz/norm_ppr.xlsx")
    data = pd.DataFrame(raw_data)
    norm_epsc = pd.pivot_table(data, values="norm_PPR",
                               index=["group", "cell"], columns="time").T
    norm_epsc.to_excel(writer)
    writer.save()

    return [info_sheet, raw_full_data, raw_full_data_min, norm_full_data, norm_full_data_min]

def personal_summary(data, good_cells=None):
    """
    just a simple little script to create a nice file for me to compare pre and
    post, both normalized and raw. EzPz baby
    """
    output = pd.DataFrame()
    if good_cells:
        cells = good_cells
    else:
        cells = data.cell.unique()
    for cell in cells:
        out = {}
        d = data.loc[data.cell == cell]
        out["cell"] = cell
        out["group"] = d.group.iloc[0]
        out["pre_epsc1"] = d.loc[(d.time > -6) & (d.time < 0),
                                 "EPSC1_amplitude"].mean()
        out["pre_epsc2"] = d.loc[(d.time > -6) & (d.time < 0),
                                 "EPSC2_amplitude"].mean()
        out["pre_ppr"] = d.loc[(d.time > -6) & (d.time < 0), "PPR"].mean()
        out["pre_norm_epsc1"] = d.loc[(d.time > -6) & (d.time < 0),
                                      "norm_EPSC1_amplitude"].mean()
        out["pre_norm_epsc2"] = d.loc[(d.time > -6) & (d.time < 0),
                                      "norm_EPSC2_amplitude"].mean()
        out["pre_norm_ppr"] = d.loc[(d.time > -6) & (d.time < 0),
                                    "norm_PPR"].mean()
        out["post_epsc1"] = d.loc[(d.time > 24) & (d.time < 30),
                                  "EPSC1_amplitude"].mean()
        out["post_epsc2"] = d.loc[(d.time > 24) & (d.time < 30),
                                  "EPSC2_amplitude"].mean()
        out["post_ppr"] = d.loc[(d.time > 24) & (d.time < 30), "PPR"].mean()
        out["post_norm_epsc1"] = d.loc[(d.time > 24) & (d.time < 30),
                                       "norm_EPSC1_amplitude"].mean()
        out["post_norm_epsc2"] = d.loc[(d.time > 24) & (d.time < 30),
                                       "norm_EPSC2_amplitude"].mean()
        out["post_norm_ppr"] = d.loc[(d.time > 24) & (d.time < 30),
                                     "norm_PPR"].mean()
        output = output.append(out, ignore_index=True)
    vars = [
        "cell",
        "group",
        "pre_epsc1",
        "post_epsc1",
        "pre_epsc2",
        "post_epsc2",
        "pre_ppr",
        "post_ppr",
        "pre_norm_epsc1",
        "post_norm_epsc1",
        "pre_norm_epsc2",
        "post_norm_epsc2",
        "pre_norm_ppr",
        "post_norm_ppr",
    ]
    output = output[vars]
    return output

def summary_graph(data, good_cells=None):
    fig, ax = plt.subplots(1,2,figsize=(10,3))
    fig.set_facecolor("w")
    textcol = "#282828"
    sns.set_style('whitegrid')
    c_strip = ["red", "black"]
    for a in ax.flatten():
        a.set_facecolor("w")
        a.spines["right"].set_visible(False)
        a.spines["top"].set_visible(False)
        a.spines["left"].set_color("#282828")
        a.spines["bottom"].set_color("#282828")

    mpl.rcParams['text.color'] = textcol
    mpl.rcParams['axes.labelcolor'] = textcol
    mpl.rcParams['xtick.color'] = textcol
    mpl.rcParams['ytick.color'] = textcol

    baseline_neg = data.loc[(data.group == "neg") & (data.time < 0)]
    baseline_pos = data.loc[(data.group == "pos") & (data.time < 0)]
    post_neg = data.loc[(data.group == "neg") & (data.time > 0)]
    post_pos = data.loc[(data.group == "pos") & (data.time > 0)]
        
    baseline_pos_mean = baseline_pos.groupby("time").mean().reset_index()
    baseline_pos_std = baseline_pos.groupby("time").sem().reset_index()
    ax[0].plot(baseline_pos_mean.time+.5, baseline_pos_mean.norm_EPSC1_amplitude,
            color=c_strip[1], lw=1,label='EAAT4+')
    ax[0].scatter(baseline_pos_mean.time+.5, baseline_pos_mean.norm_EPSC1_amplitude,
               color=c_strip[1], s=15, edgecolor="#282828")
    ax[0].errorbar(baseline_pos_mean.time+.5, baseline_pos_mean.norm_EPSC1_amplitude,
                yerr=baseline_pos_std.norm_EPSC1_amplitude, fmt="none",
                color=c_strip[1], lw=.5)
    
    baseline_neg_mean = baseline_neg.groupby("time").mean().reset_index()
    baseline_neg_std = baseline_neg.groupby("time").sem().reset_index()
    ax[0].plot(baseline_neg_mean.time, baseline_neg_mean.norm_EPSC1_amplitude,
            color=c_strip[0], lw=1,label='EAAT4-')
    ax[0].scatter(baseline_neg_mean.time, baseline_neg_mean.norm_EPSC1_amplitude,
               color=c_strip[0], s=15, edgecolor="#282828")
    ax[0].errorbar(baseline_neg_mean.time, baseline_neg_mean.norm_EPSC1_amplitude,
                yerr=baseline_neg_std.norm_EPSC1_amplitude, fmt="none",
                color=c_strip[0], lw=.5)

    post_neg_mean = post_neg.groupby("time").mean().reset_index()
    post_neg_std = post_neg.groupby("time").std().reset_index()
    ax[0].plot(post_neg_mean.time, post_neg_mean.norm_EPSC1_amplitude,
            color=c_strip[0], lw=1)
    ax[0].scatter(post_neg_mean.time, post_neg_mean.norm_EPSC1_amplitude,
               color=c_strip[0], s=15, edgecolor="#282828")
    ax[0].errorbar(post_neg_mean.time, post_neg_mean.norm_EPSC1_amplitude,
                yerr=post_neg_std.norm_EPSC1_amplitude, fmt="none",
                color=c_strip[0], lw=.5)
    
    post_pos_mean = post_pos.groupby("time").mean().reset_index()
    post_pos_std = post_pos.groupby("time").std().reset_index()
    ax[0].plot(post_pos_mean.time+.5, post_pos_mean.norm_EPSC1_amplitude,
            color=c_strip[1], lw=1)
    ax[0].scatter(post_pos_mean.time+.5, post_pos_mean.norm_EPSC1_amplitude,
               color=c_strip[1], s=15, edgecolor="#282828")
    ax[0].errorbar(post_pos_mean.time+.5, post_pos_mean.norm_EPSC1_amplitude,
                yerr=post_pos_std.norm_EPSC1_amplitude, fmt="none",
                color=c_strip[1], lw=.5)
    ax[0].set_ylim(0,2)
    
    ax[1].plot(baseline_neg_mean.time, baseline_neg_mean.norm_PPR,
            color=c_strip[0], lw=1)
    ax[1].scatter(baseline_neg_mean.time, baseline_neg_mean.norm_PPR,
               color=c_strip[0], s=15, edgecolor="#282828")
    ax[1].errorbar(baseline_neg_mean.time, baseline_neg_mean.norm_PPR,
                yerr=baseline_neg_std.norm_PPR, fmt="none",
                color=c_strip[0], lw=.5)
    
    post_pos_mean = post_pos.groupby("time").mean().reset_index()
    post_pos_std = post_pos.groupby("time").std().reset_index()
    ax[1].plot(post_pos_mean.time+.5, post_pos_mean.norm_PPR,
            color=c_strip[1], lw=1)
    ax[1].scatter(post_pos_mean.time+.5, post_pos_mean.norm_PPR,
               color=c_strip[1], s=15, edgecolor="#282828", label="EAAT4+")
    ax[1].errorbar(post_pos_mean.time+.5, post_pos_mean.norm_PPR,
                yerr=post_pos_std.norm_PPR, fmt="none",
                color=c_strip[1], lw=.5)
    
    post_neg_mean = post_neg.groupby("time").mean().reset_index()
    post_neg_std = post_neg.groupby("time").std().reset_index()
    ax[1].plot(post_neg_mean.time, post_neg_mean.norm_PPR,
            color=c_strip[0], lw=1)
    ax[1].scatter(post_neg_mean.time, post_neg_mean.norm_PPR,
               color=c_strip[0], s=15, edgecolor="#282828", label="EAAT4-")
    ax[1].errorbar(post_neg_mean.time, post_neg_mean.norm_PPR,
                yerr=post_neg_std.norm_PPR, fmt="none",
                color=c_strip[0], lw=.5)

    
    baseline_pos_mean = baseline_pos.groupby("time").mean().reset_index()
    baseline_pos_std = baseline_pos.groupby("time").sem().reset_index()
    ax[1].plot(baseline_pos_mean.time+.5, baseline_pos_mean.norm_PPR,
            color=c_strip[1], lw=1)
    ax[1].scatter(baseline_pos_mean.time+.5, baseline_pos_mean.norm_PPR,
               color=c_strip[1], s=15, edgecolor="#282828")
    ax[1].errorbar(baseline_pos_mean.time+.5, baseline_pos_mean.norm_PPR,
                yerr=baseline_pos_std.norm_PPR, fmt="none",
                color=c_strip[1], lw=.5)

    
    ax[1].set_ylim(0,2)
    ax[0].legend(loc='upper left')
    ax[1].set_visible(False)
    for a in ax.flatten():
        a.set_xlim(-5.5,35.5)
        a.set_ylim(.5,1.65)
        a.set_xlabel("time (min)")
    ax[0].set_ylabel("Normalized EPSC1 amplitude")
    ax[1].set_ylabel("normalized PPR")
    ax[0].set_title('A',loc='left', fontweight ='bold', fontsize = '16')
    fig.tight_layout()
    return fig, ax

def complete_plot(data):
    plot_data  = data.groupby(['group','cell','pre_post']).median().reset_index()
    temp_data  = plot_data.groupby('pre_post')
    pre_data   = temp_data.get_group('pre')
    post_data  = temp_data.get_group('post')
    
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(2,3)
    
    vars = ['norm_EPSC1_amplitude','norm_EPSC2_amplitude','PPR']
    for i,v in enumerate(vars):
        sns.boxplot(x='group',y=v,data= pre_data,ax=ax[0,i], fliersize=0, order=['pos','neg'])
        sns.stripplot(x='group',y=v,data= pre_data,ax=ax[0,i], palette=['black'], jitter=0.1, s=3, order=['pos','neg'])
        
        sns.boxplot(x='group',y=v,data= post_data,ax=ax[1,i], fliersize=0, order=['pos','neg'])
        sns.stripplot(x='group',y=v,data= post_data,ax=ax[1,i], palette=['black'], jitter=0.1, s=3, order=['pos','neg'])
        
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.set_xlabel('')
        
        
    ax[0,0].set_ylabel('Normalized EPSC1 pre')
    ax[0,1].set_ylabel('Normalized EPSC2 pre')
    ax[1,0].set_ylabel('Normalized EPSC1 post')
    ax[1,1].set_ylabel('Normalized EPSC2 post')
    ax[0,2].set_ylabel('PPR pre')
    ax[1,2].set_ylabel('PPR post')
    ax[0,0].set_xticklabels([' ', ' '])
    ax[0,1].set_xticklabels([' ', ' '])
    ax[0,2].set_xticklabels([' ', ' '])
    ax[1,0].set_xticklabels(['EAAT4+','EAAT4-'])
    ax[1,1].set_xticklabels(['EAAT4+','EAAT4-'])
    ax[1,2].set_xticklabels(['EAAT4+','EAAT4-'])
    
    ax[0,0].set_title('A', loc='left')
    ax[0,1].set_title('C', loc='left')
    ax[0,2].set_title('E', loc='left')
    ax[1,0].set_title('B', loc='left')
    ax[1,1].set_title('D', loc='left')
    ax[1,2].set_title('F', loc='left')
    
    ax[0,0].set_yticks([0.9,0.95,1.0,1.05,1.1])
    ax[0,1].set_yticks([0.9,0.95,1.0,1.05,1.1])
    
    ax[1,0].set_yticks([0.75,1.00,1.25,1.50,1.75,2.00])
    ax[1,1].set_yticks([0.75,1.00,1.25,1.50,1.75,2.00])
    
    ax[0,2].set_yticks([1.0,1.2,1.4,1.6,1.8,2.0])
    ax[1,2].set_yticks([1.0,1.2,1.4,1.6,1.8])
    
    fig.tight_layout()
    return fig, ax

def ltp_stats(data):
    plot_data  = data.groupby(['group','cell','pre_post']).median().reset_index()
    temp_data  = plot_data.groupby('pre_post')
    pre_data   = temp_data.get_group('pre')
    post_data  = temp_data.get_group('post')
    
    vars = ['norm_EPSC1_amplitude','norm_EPSC2_amplitude','PPR']
    for v in vars:
        print(v, tstd(pre_data.loc[pre_data.group == "pos", v]), 'Z+ standard deviation')
        print(v, tstd(pre_data.loc[pre_data.group == "neg", v]), 'Z- standard deviation')
        print(v, levene(pre_data.loc[pre_data.group == "pos", v], pre_data.loc[pre_data.group == "neg", v]))
        print(v, ttest_ind(pre_data.loc[pre_data.group == "pos", v], pre_data.loc[pre_data.group == "neg", v]), "\n")
        
    for v in vars:
        print(v, tstd(post_data.loc[post_data.group == "pos", v]), 'Z+ standard deviation')
        print(v, tstd(post_data.loc[post_data.group == "neg", v]), 'Z- standard deviation')
        print(v, levene(post_data.loc[post_data.group == "pos", v], post_data.loc[post_data.group == "neg", v]))
        print(v, ttest_ind(post_data.loc[post_data.group == "pos", v], post_data.loc[post_data.group == "neg", v]), "\n")
        
    return
    
if __name__ == "__main__":
    base_folder = "C:/Users/Gebruiker/Desktop/EAAT4_IPSC_EPSC/LTP_apr"
    data_loc = os.path.join(base_folder, "ABF/")
    meta_loc = os.path.join(base_folder, "meta/ltpmeta_select.xlsx")
    output_loc = os.path.join(base_folder, "output/")
    # output = find_meta(meta_loc, data_loc)


