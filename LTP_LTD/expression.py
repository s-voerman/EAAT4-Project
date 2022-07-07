import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.interpolate import make_interp_spline

def plcb4_correlation():
    os.chdir(data_loc)
    eaat4 = pd.read_csv('EAAT4.csv')
    aldoc = pd.read_csv('Aldoc.csv')
    plcb4 = pd.read_csv('PLCB4.csv')
    
    #normalization
    eaat4['norm_gray_value'] = eaat4.Gray_Value / np.mean(eaat4.Gray_Value)
    aldoc['norm_gray_value'] = aldoc.Gray_Value / np.mean(aldoc.Gray_Value)
    plcb4['norm_gray_value'] = plcb4.Gray_Value / np.mean(plcb4.Gray_Value)
     
    cm=1/2.54
    fig, ax = plt.subplots(3,1, figsize=(9.15*cm, 5 *cm))
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.axhline(y=1, lw=0.5, linestyle='--', color='gray')
    sns.lineplot(x='Distance_(microns)', y='norm_gray_value', data=eaat4, ax=ax[0], color='#25A73B')
    sns.lineplot(x='Distance_(microns)', y='norm_gray_value', data=plcb4, ax=ax[1], color='#D42F2F')
    sns.lineplot(x='Distance_(microns)', y='norm_gray_value', data=aldoc, ax=ax[2], color='#2596be')
    
    for a in ax.flatten():
        a.set_ylabel(' ')
    ax[0].spines['bottom'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[0].axes.get_xaxis().set_visible(False)
    ax[1].axes.get_xaxis().set_visible(False)
    ax[0].set_ylim(0,2.0)
    ax[1].set_ylim(0.5,1.5)
    ax[2].set_ylim(0.5,1.5)
    ax[0].set_yticks([0,1,2])
    ax[1].set_yticks([0.5,1,1.5])
    ax[2].set_yticks([0.5,1,1.5])
    ax[0].set_yticklabels([0,100,200])
    ax[1].set_yticklabels([50,100,150])
    ax[2].set_yticklabels([50,100,150])
    ax[1].set_ylabel('Relative expression to the mean')
    ax[0].set_xlabel('')
    ax[1].set_xlabel('')
    ax[2].set_xlabel('Distance (microns)')
    
    #plt.subplots_adjust(wspace=0, hspace=0)
    #fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'expression.pdf'))
    return 
    
if __name__ == "__main__":
    base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
    data_loc = os.path.join(base_folder, "Data")
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
