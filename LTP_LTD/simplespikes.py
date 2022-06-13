import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from neo.rawio import AxonRawIO
import os

def read_abf(file_loc):
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header["signal_channels"][0][2])
    data = a.get_analogsignal_chunk() * a.header["signal_channels"][0][5]
    return data[:, 0], sr

def example_trace():
    os.chdir(data_loc)
    data, sr = read_abf('cclamp_2.ABF')
    
    fig, ax = plt.subplots()
    
    ax.plot(data[0:10000], color='black')
    
    fig.savefig(os.path.join(output_loc, "simplespikes.pdf"), dpi = 300)
    return fig, ax

if __name__ == "__main__":
    base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
    data_loc = os.path.join(base_folder, "Data/simplespikes")
    output_loc = os.path.join(base_folder, "Analysis_tools/output/")
    plt.style.use('seaborn-paper')
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
