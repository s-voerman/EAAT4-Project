import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from neo.rawio import AxonRawIO

base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
ltd_loc = os.path.join(base_folder, "Data/LTD/Data")
ltp_loc = os.path.join(base_folder, "Data/LTP/Combined/abf")
ltp_loc_pt = os.path.join(base_folder, "Data/LTP/PT")
mepsc_loc = os.path.join(base_folder, "Data/mEPSCs/Combined/epsc")
mipsc_loc = os.path.join(base_folder, "Data/mIPSCs/Combined/ipsc")
meta_loc = os.path.join(base_folder, "Metadata/ltdmeta.xlsx")
output_loc = os.path.join(base_folder, "Analysis_tools/output/")

def read_abf(file_loc):
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header["signal_channels"][0][2])
    data = a.get_analogsignal_chunk() * a.header["signal_channels"][0][5]
    return data[:, 0], sr

def example_trace():
    os.chdir(mipsc_loc)
    
    #LTD
    # neg_data_pre, sr = read_abf('stijn_120122 1-58.ABF')
    # neg_data_post, sr = read_abf('stijn_120122 1-60.ABF')
    # pos_data_pre, sr = read_abf('stijn_190122 3-101.ABF')
    # pos_data_post, sr = read_abf('stijn_190122 3-103.ABF')
    
    #LTP
    # pos_data_pre, sr = read_abf('stijn_030720 1-18.ABF')
    # pos_data_post, sr = read_abf('stijn_030720 1-20.ABF')
    # neg_data_pre, sr = read_abf('stijn_270820 1-9.ABF')
    # neg_data_post, sr = read_abf('stijn_270820 1-11.ABF')
    
    #LTP_PT
    # pos_data_pre, sr = read_abf('stijn_010322 2-36.ABF')
    # pos_data_post, sr = read_abf('stijn_010322 2-38.ABF')
    # neg_data_pre, sr = read_abf('stijn_240222 1-31.ABF')
    # neg_data_post, sr = read_abf('stijn_240222 1-33.ABF')
    
    #Flocculus
    # pre_data,sr = read_abf('stijn_060422 1-32.ABF')
    # post_data,sr  = read_abf('stijn_060422 1-34.ABF')
    
    #mEPSC
    # pos_data, sr = read_abf('stijn_260121 3-14.ABF') 
    # neg_data, sr = read_abf('stijn_120320 1-3.ABF')
    
    # #mIPSC
    pos_data, sr = read_abf('laura_271020 5-27.ABF') 
    neg_data, sr = read_abf('laura_100920 5-19.ABF')
    
    #Plot
    cm = 1/2.54
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rcParams['font.family'] = 'sans-serif'
    fonts = {"font.size":7, "axes.labelsize":7, "ytick.labelsize":7, "xtick.labelsize":7}
    plt.rcParams.update(fonts)
    fig, ax = plt.subplots(1,2, dpi=300)
    
    ax[0].plot(pos_data[0:50000])
    ax[1].plot(neg_data[0:50000])
 

    
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'mipsc_example_trace.png'), dpi=300)

    return fig, ax


