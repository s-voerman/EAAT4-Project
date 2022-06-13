import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from scipy.stats import levene, ttest_ind
from scipy.stats import tstd

plt.style.use('seaborn-paper')
sns.set_palette('Set1',8,0.7)
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['xtick.minor.width'] = 0.3
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['ytick.minor.width'] = 0.3
fonts = {"font.size":7, "axes.labelsize":7, "ytick.labelsize":7, "xtick.labelsize":7}
plt.rcParams.update(fonts)

#PLOTTING

def plot_sealdata():
    #Combined data
    plot_data = pd.read_excel('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/seal_data_combined.xlsx')
    plot_data.loc[:,'Ratio (Rs/Ri)'] = (plot_data.Rs / plot_data.Ri) *100
    
    cm=1/2.54
    fig, ax     = plt.subplots(1,3, figsize=(9.15*cm,3))    
    vars        = ['Rm','Cm','Tau']
    for i,v in enumerate(vars):
        sns.boxplot(data=plot_data, x='group', y=v, ax=ax[i], fliersize=0, order=['pos','neg'], linewidth=0.5)
        sns.stripplot(data=plot_data, x='group', y=v, ax=ax[i], palette=["#000000","#FFFFFF"], edgecolor='black', linewidth=0.5
                      ,order=['pos','neg'], size=2, marker='^')
    
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.set_xticklabels(['+','-'])
        a.set_xlabel('EAAT4')
    ax[0].set_ylim(0,500)
    ax[0].set_ylabel('Membrane resistance (MΩ)')
    ax[1].set_ylim(0,2000)
    ax[1].set_ylabel('Membrane capacitance (pF)')
    ax[2].set_ylim(0,500)
    ax[2].set_ylabel('Tau (ms)')

    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "seal_plots.pdf"))
    
    return fig, ax
        
def lobule_plot():
    plot_data = pd.read_excel('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/seal_data_combined.xlsx')
    
    cm=1/2.54
    fig, ax = plt.subplots(dpi=300, figsize=(11.4 *cm,8.3*cm))
    sns.countplot(data=plot_data, x='lobule', ax=ax, hue='group', hue_order=['pos','neg'], saturation=.5,
                  order=['II','III','IV-V','VI','VII','VIII','IX','Crus1','Paramedian','Simple','?'])
    
    ax.set_xticklabels(['II','III','IV-V','VI','VII','VIII','IX', 'Crus1','PM','SM','Unid'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0,30)
    ax.set_xlabel(' ')
    ax.set_ylabel('n')
    fig.tight_layout()
    fig.savefig('celllocs.png',facecolor = 'white', edgecolor='white',dpi=1000)
    
    return fig, ax

#STATS

def seal_stats():
    excel_df = pd.read_excel('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/seal_data_combined.xlsx')
    
    vars        = ['Rm','Cm','Tau']
    for v in vars:
        print(v, tstd(excel_df.loc[excel_df.group == "pos", v]), 'Z+ standard deviation')
        print(v, tstd(excel_df.loc[excel_df.group == "neg", v]), 'Z- standard deviation')
        print(v, levene(excel_df.loc[excel_df.group == "pos", v], excel_df.loc[excel_df.group == "neg", v]))
        print(v, ttest_ind(excel_df.loc[excel_df.group == "pos", v], excel_df.loc[excel_df.group == "neg", v]), "\n")
 
    
def domain_stats():
    excel_df = pd.read_excel('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/seal_data_combined.xlsx')
    
    vars = ['Rs','Rm','cslow','holding','Ratio (Rs/Ri)','Cm']
    for v in vars:
        print(v, tstd(excel_df.loc[excel_df.domain == 0, v]), 'AZ standard deviation')
        print(v, tstd(excel_df.loc[excel_df.domain == 2, v]), 'PZ standard deviation')
        print(v, levene(excel_df.loc[excel_df.domain == 0, v], excel_df.loc[excel_df.domain == 2, v]))
        print(v, ttest_ind(excel_df.loc[excel_df.domain == 0, v], excel_df.loc[excel_df.domain == 2, v]), "\n")

def zebrin_stats():
    excel_df = pd.read_excel('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/seal_data_combined.xlsx')
    AZ_df = excel_df.loc[excel_df.domain == 0]
    CZ_df = excel_df.loc[excel_df.domain == 1]
    PZ_df = excel_df.loc[excel_df.domain == 2]
    
    vars = ['Ri','Rs','Rm','cslow','holding','Ratio (Rs/Ri)','Cm']
    for v in vars:
        
        print(v, tstd(AZ_df.loc[AZ_df.group == "pos", v]), 'AZ Z+ standard deviation')
        print(v, tstd(AZ_df.loc[AZ_df.group == "neg", v]), 'AZ Z- standard deviation')
        print(v, levene(AZ_df.loc[AZ_df.group == "pos", v], AZ_df.loc[AZ_df.group == "neg", v]))
        print(v, ttest_ind(AZ_df.loc[AZ_df.group == "pos", v], AZ_df.loc[AZ_df.group == "neg", v]), "\n")
        
        print(v, tstd(CZ_df.loc[CZ_df.group == "pos", v]), 'CZ Z+ standard deviation')
        print(v, tstd(CZ_df.loc[CZ_df.group == "neg", v]), 'CZ Z- standard deviation')
        print(v, levene(CZ_df.loc[CZ_df.group == "pos", v], CZ_df.loc[CZ_df.group == "neg", v]))
        print(v, ttest_ind(CZ_df.loc[CZ_df.group == "pos", v], CZ_df.loc[CZ_df.group == "neg", v]), "\n")
        
        print(v, tstd(PZ_df.loc[PZ_df.group == "pos", v]), 'PZ Z+ standard deviation')
        print(v, tstd(PZ_df.loc[PZ_df.group == "neg", v]), 'PZ Z- standard deviation')
        print(v, levene(PZ_df.loc[PZ_df.group == "pos", v], PZ_df.loc[PZ_df.group == "neg", v]))
        print(v, ttest_ind(PZ_df.loc[PZ_df.group == "pos", v], PZ_df.loc[PZ_df.group == "neg", v]), "\n")
    
    
if __name__ == "__main__":
    base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
    data_loc = os.path.join(base_folder, "Data/LTD/Data")
    meta_loc = os.path.join(base_folder, "Metadata/ltdmeta_select.xlsx")
    output_loc = os.path.join(base_folder, "Analysis_tools/output/")
    plt.style.use('seaborn-paper')
    color = ["#5BE12D","#D24646"] 
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

    
    
    
    
    
    
    