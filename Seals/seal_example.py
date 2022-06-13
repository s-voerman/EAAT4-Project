import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from neo.rawio import AxonRawIO

def read_abf(file_loc):
    #os.chdir('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/LTD')
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header["signal_channels"][0][2])
    data = a.get_analogsignal_chunk() * a.header["signal_channels"][0][5]
    
    return data[:, 0], sr

def seal_plot():
    os.chdir('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/mEPSCs/Combined/seals/abf')
    data,sr = read_abf('stijn_120320 1-2.ABF')
    
    cm=1/2.54
    fig, ax = plt.subplots(dpi=300, figsize=(9.15*cm,4))
    ax.plot(data[int(0.65*sr):int(0.85*sr)], color='black')
    
    os.chdir(output_loc)
    fig.savefig('seal_example.png', dpi=1000)
    return fig, ax

def ltd_plot():
    pre_data, sr = read_abf('stijn_121121 2-83.ABF')
    post_data, sr = read_abf('stijn_121121 2-85.ABF')
    
    fig, ax = plt.subplots(dpi=300, figsize = (4,4))
    ax.plot(pre_data[0:20000], label='pre')
    ax.plot(post_data[100000:120000]+20, label='post')
    plt.legend()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    fig.savefig('LTD_example', facecolor='white')
    return fig, ax 

def indplot():
    data,sr = read_abf('stijn_091121 2-56.ABF')
    fig, ax = plt.subplots(dpi = 300)
    
    ax.plot(data[0:20000])
    
    fig.savefig('LTD_induction.png', facecolor='white')
    return fig, ax

def epsc_plot():
    os.chdir('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/mEPSCs/By_date/270220')
    data, sr = read_abf('stijn_270220 4-24.ABF')
    
    fig, ax = plt.subplots(1,1)
    ax.plot(data[3850000:3900000], color='black')
    ax.set_ylim(-200,0)
    
    return fig, ax

if __name__ == "__main__":
    base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
    data_loc = os.path.join(base_folder, "Data/LTD/Data")
    meta_loc = os.path.join(base_folder, "Metadata/ltdmeta.xlsx")
    output_loc = os.path.join(base_folder, "Analysis_tools/output/")