import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from neo.rawio import AxonRawIO
from numpy import trapz
from scipy.stats import levene, ttest_ind
from scipy.stats import tstd, sem, shapiro, mannwhitneyu

#Essential variables with directories of the data, metafiles, and output location, change this to whereever you save these items
base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
data_loc = os.path.join(base_folder, "Data/LTD/Data")
meta_loc = os.path.join(base_folder, "Metadata/ltdmeta_select_strict.xlsx")
output_loc = os.path.join(base_folder, "Analysis_tools/output/")

#Simple function that reads abf files and returns data as a tuple containing data (list) and sampling rate (int)
def read_abf(file_loc):
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header["signal_channels"][0][2])
    data = a.get_analogsignal_chunk() * a.header["signal_channels"][0][5]
    return data[:, 0], sr

def meta():
    meta = pd.read_excel(meta_loc)
    for i, x in meta.iterrows():
        fn = x["name"] + " " + str(x.cell) + "-" + str(x.ind)
        fn = fn + ".ABF"
        meta.loc[i, "filename"] = fn
    return meta

#Analyzes complex spikes based on the list of files in meta
def complex_spikes():
    cs_meta = meta() #reads meta
    os.chdir(data_loc) #changes working directory to location of data
    for i, x in cs_meta.iterrows():
        data, sr = read_abf(x.filename) #gets the data
        window = np.arange(0, len(data), sr) #determines the window of 1 sweep
        control_df = pd.DataFrame() #output frame
        Area = [] #list for temporary storage
        for start in window:
            baseline = np.min(data[int(0.3 * sr)+start:int(0.4 * sr)+start]) #determines the baseline holding voltage
            cs_data = data[int(0.20 * sr)+start:int(0.24 * sr)+start] #gets the data of the complex spike
            cs_data = cs_data- baseline #normalizes everything to be at V = 0
            x_data = np.arange(0,40,0.02) #x_data based on the length of the complex spike (in this case 40 ms)
            area = trapz(cs_data, x_data) #determines the area of the complex spike from cs_data and x_data
            Area.append(area)
            # plt.plot(x_data, cs_data, color='black')
            # plt.fill_between(x_data, cs_data,alpha=0.25)
            # plt.savefig(os.path.join(output_loc, "summary_" + str(x.name) + ".png"))
            # plt.close()
        control_df['area'] = pd.Series(Area)
        cs_meta.loc[i,'area'] = np.mean(Area)
    return cs_meta

#plots a single cell's complex spike
def cs_plot():
    os.chdir(data_loc)
    data, sr = read_abf('stijn_141221 1-157.ABF')
    window = np.arange(0, len(data), sr)
    Area = []
    for start in window:
        baseline = np.min(data[int(0.3 * sr)+start:int(0.4 * sr)+start])
        cs_data = data[int(0.05 * sr)+start:int(0.30 * sr)+start]
    plt.plot(cs_data)
    plt.savefig('cs_test.png', dpi=1000, color='black')
    return Area

#plots the area of all cells based on pos/neg
def plot_area(cs_meta):
    cm= 1/2.54
    fig, ax = plt.subplots(figsize=(3.4*cm,3))
    sns.boxplot(x='group',y='area', data=cs_meta, fliersize=0, ax=ax, linewidth=0.5)
    sns.swarmplot(x='group',y='area', data=cs_meta, ax=ax, palette=["#000000","#FFFFFF"], edgecolor='black', linewidth=0.5, size=3, 
                  marker = '^')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0,500)
    ax.set_xticklabels(['+','-'])
    ax.set_xlabel('EAAT4')
    ax.set_ylabel('Area (mV x ms)')
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "complex_spikedata.pdf"))
    return fig, ax

#plots a single cs, for figures
def stim(cs_meta):
    os.chdir(data_loc)
    data, sr = read_abf('stijn_141221 1-157.ABF')
    fig, ax = plt.subplots(dpi=300)
    stim_data = data[int(0.09*sr):int(0.26*sr)]
    ax.plot(stim_data,color='black')
    return

#A really bloated function that does statistics for complex spike area. Not my best work.
def stats(cs_meta):
    meta = {}
    vars = ['area']
    pos_data = cs_meta.loc[cs_meta.group == 'pos']
    neg_data = cs_meta.loc[cs_meta.group == 'neg']
    for v in vars:
        neg_mean = np.mean(neg_data.loc[:,v])
        pos_mean = np.mean(pos_data.loc[:,v])
        neg_sem = sem(neg_data.loc[:,v])
        pos_sem = sem(pos_data.loc[:,v])
        neg_sd = tstd(neg_data.loc[:,v])
        pos_sd = tstd(pos_data.loc[:,v])
        pos_sh = shapiro(pos_data.loc[:,v])
        neg_sh = shapiro(neg_data.loc[:,v])
        meta['pos_mean'] = pos_mean
        meta['neg_mean'] = neg_mean
        meta['pos_sem'] = pos_sem
        meta['neg_sem'] = neg_sem
        meta['pos_sd'] = pos_sd
        meta['neg_sd'] = neg_sd
        meta['pos_Sh-W'] = pos_sh[1]
        meta['neg_Sh-W'] = neg_sh[1]
        print(v, tstd(cs_meta.loc[cs_meta.group == "pos", v]), 'Z+ standard deviation')
        print(v, tstd(cs_meta.loc[cs_meta.group == "neg", v]), 'Z- standard deviation')
        print(v, levene(cs_meta.loc[cs_meta.group == "pos", v], cs_meta.loc[cs_meta.group == "neg", v]))
        print(v, ttest_ind(cs_meta.loc[cs_meta.group == "pos", v], cs_meta.loc[cs_meta.group == "neg", v]), "\n")
        print(v, mannwhitneyu(cs_meta.loc[cs_meta.group == "pos", v], cs_meta.loc[cs_meta.group == "neg", v]), "\n")
    stat_df = pd.DataFrame([meta])
    stat_df.to_excel(os.path.join(output_loc, 'LTD ' + str(v) + '.xlsx'))
    return stat_df

#Variable definitions for plotting
if __name__ == "__main__":
    plt.style.use('seaborn-paper')
    color = ["#5BE12D","#EC008C"] 
    sns.set_palette(sns.color_palette(color), desat=1)
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['xtick.minor.width'] = 0.3
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['ytick.minor.width'] = 0.3
    fonts = {"font.size":7, "axes.labelsize":7, "ytick.labelsize":7, "xtick.labelsize":7}
    plt.rcParams.update(fonts)
