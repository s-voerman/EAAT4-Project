import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import math
import os

def plot(data):
    pre_neg = data.loc[(data.group == "neg") & (data.time < 0) & (data.time > -11)]
    pre_pos = data.loc[(data.group == "pos") & (data.time < 0) & (data.time > -11)]
    post_neg = data.loc[(data.group == "neg") & (data.time >0) & (data.time < 36)]
    post_pos = data.loc[(data.group == "pos") & (data.time >0) & (data.time < 36)]
    
    pre_neg_std = pre_neg.groupby('time').std().reset_index()
    pre_pos_std = pre_pos.groupby('time').std().reset_index()
    post_neg_std = post_neg.groupby('time').std().reset_index()
    post_pos_std = post_pos.groupby('time').std().reset_index()
    n_pos_cells = 13
    n_neg_cells = 13
    yerror_1 = pre_neg_std.norm_EPSC1_amplitude / math.sqrt(n_neg_cells)
    yerror_2 = pre_pos_std.norm_EPSC1_amplitude / math.sqrt(n_pos_cells)
    yerror_3 = post_neg_std.norm_EPSC1_amplitude / math.sqrt(n_neg_cells)
    yerror_4 = post_pos_std.norm_EPSC1_amplitude / math.sqrt(n_pos_cells)
    
    pre_neg_mean = pre_neg.groupby('time').mean().reset_index()
    pre_pos_mean = pre_pos.groupby('time').mean().reset_index()
    post_neg_mean = post_neg.groupby('time').mean().reset_index()
    post_pos_mean = post_pos.groupby('time').mean().reset_index()

    cm = 1/2.54
    fig, ax = plt.subplots(figsize=((0.7*18.3)*cm,3), dpi=300)
    ax.errorbar('time','norm_EPSC1_amplitude', data=pre_neg_mean, yerr= yerror_1, fmt=' ', elinewidth=0.5, capsize=2.0, capthick=0.5,)
    ax.errorbar('time','norm_EPSC1_amplitude', data=pre_pos_mean, yerr= yerror_2, fmt=' ', elinewidth=0.5, capsize=2.0, capthick=0.5)
    ax.errorbar('time','norm_EPSC1_amplitude', data=post_neg_mean, yerr= yerror_3, fmt=' ', elinewidth=0.5, capsize=2.0, capthick=0.5)
    ax.errorbar('time','norm_EPSC1_amplitude', data=post_pos_mean, yerr= yerror_4, fmt=' ', elinewidth=0.5, capsize=2.0, capthick=0.5)
    ax.plot('time','norm_EPSC1_amplitude', marker='o',data=pre_neg_mean, mec='black', mew=0.5)
    ax.plot('time','norm_EPSC1_amplitude', marker='^',data=pre_pos_mean, mec='black', mew=0.5)
    ax.plot('time','norm_EPSC1_amplitude', marker='o',data=post_neg_mean, mec='black', mew=0.5)
    ax.plot('time','norm_EPSC1_amplitude', marker='^',data=post_pos_mean, mec='black', mew=0.5)
    

    for a in ax.lines:
        a.set_linewidth(0.5)
    return fig, ax

if __name__ == "__main__":
    base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
    data_loc = os.path.join(base_folder, "Data/LTD/Data")
    meta_loc = os.path.join(base_folder, "Metadata/ltdmeta_select.xlsx")
    output_loc = os.path.join(base_folder, "Analysis_tools/output/")
    plt.style.use('seaborn-paper')
    #color = ["#FFFFFF","#808080"] #white/grey
    color = ["#5BE12D","#D24646"] #Original green/red
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