import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
from scipy.stats import levene, ttest_ind
from scipy.stats import linregress
from scipy.stats import tstd

########### EPSCs #################################################################################################################

def EPSC_meta():
    os.chdir('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Metadata') #Change this to the folder with the relevant metafile
    meta = pd.read_excel('mepsc_metafile.xlsx')
    meta['filename'] = [" ".join([x, "-".join([str(y),str(z)])]) for x,y,z in zip(meta.name, meta.cell, meta.epsc)]
    vars = ['II','III','IV-V'] #AZ
    for i,v in enumerate(vars):
        meta.loc[meta['lobule'] == v,'domain'] = 0
    vars = ['VI','VII'] #CZ
    for i,v in enumerate(vars):
        meta.loc[meta['lobule'] == v,'domain'] = 1
    vars = ['VIII','IX'] #PZ
    for i,v in enumerate(vars):
        meta.loc[meta['lobule'] == v,'domain'] = 2
    return meta

def EPSC_reader():
    meta = EPSC_meta()
    os.chdir('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/mEPSCs/Combined/excel/clean') #Directory that contains EPSC !Excel! files
    output_df = pd.DataFrame()
    dirlist = os.listdir() 
    
    for filename in dirlist:
        if not filename.endswith('.xlsx'):
            continue 
        cellname = filename.split(".")[0]
        print(meta.loc[meta.filename == filename.split(".")[0], "group"])
        print(cellname)
        zebrin_id = meta.loc[meta.filename == filename.split(".")[0],"group"].item()
        lobule = meta.loc[meta.filename == filename.split(".")[0],"lobule"].item()
        confirmed = meta.loc[meta.filename == filename.split(".")[0],"confirmed"].item()
        domain = meta.loc[meta.filename == filename.split(".")[0],"domain"].item()
        df = pd.read_excel(filename)
        df = df.drop(df.columns[0], axis=1)
        Frequency = len(df) / 120
        df['Zebrin']            = zebrin_id
        df['domain']            = domain
        df['Frequency (Hz)']    = Frequency
        df['Cellname']          = filename
        df['lobule']            = lobule
        df['confirmed']         = confirmed
        df['amplitude']         = df['amplitude'] * -1 #flips all negative values to positive values
        df['rise']              = df['rise'] * 1000 #converts seconds to miliseconds
        df['tau']               = df['tau'] * 1000 #converts seconds to miliseconds
        output_df = output_df.append(df, ignore_index = True)
    return(output_df)

def exclusion_list_epsc(e):
    #creates a list of cells to exclude based on criteria
    exclusion_df = pd.DataFrame()
    
    #baseline
    baseline_cutoff = -500
    baseline_df     = e[e['baseline'] < baseline_cutoff]
    
    #frequency
    frequency_cutoff1 = 15
    frequency_cutoff2 = 0.5
    frequency_df1      = e[e['Frequency (Hz)'] > frequency_cutoff1]
    frequency_df2      = e[e['Frequency (Hz)'] < frequency_cutoff2]
    
    #rise time
    rise_cutoff     = 8 
    rise_df         = e[e['rise'] > rise_cutoff]
    
    #amplitude
    amp_cutoff      = 30 
    amp_df          = e[e['amplitude'] > amp_cutoff]
    
    #tau
    tau_cutoff      = 14
    tau_df          = e[e['tau'] > tau_cutoff]
    
    #putting it all together
    exclusion_df = exclusion_df.append(baseline_df)
    exclusion_df = exclusion_df.append(frequency_df1)
    exclusion_df = exclusion_df.append(frequency_df2)
    exclusion_df = exclusion_df.append(rise_df)
    exclusion_df = exclusion_df.append(tau_df)
    exclusion_df = exclusion_df.append(amp_df)
    exclusion_df = exclusion_df.drop_duplicates()

    return exclusion_df

def export_EPSC_data():
    output_df = EPSC_reader()
    temp_median = output_df.groupby(['Zebrin','Cellname']).median().reset_index()
    temp_mean = output_df.groupby(['Zebrin','Cellname']).mean().reset_index()
    
    os.chdir('//storage.erasmusmc.nl/m/MyDocs/051071/My Documents/Desktop/Analysis_tools/output')
    temp_median.to_excel('./EPSC_mediandata.xlsx')
    temp_mean.to_excel('./EPSC_meandata.xlsx')
    return


####### IPSCs ####################################################################################################################

def IPSC_meta():
    os.chdir('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Metadata') #Change this to the folder with the relevant metafile
    meta = pd.read_excel('mipsc_metafile.xlsx')
    meta['filename'] = [" ".join([x, "-".join([str(y),str(z)])]) for x,y,z in zip(meta.name, meta.cell, meta.ipsc)]
    
    vars = ['II','III','IV-V'] #AZ
    for i,v in enumerate(vars):
        meta.loc[meta['lobule'] == v,'domain'] = 0
    vars = ['VI','VII'] #CZ
    for i,v in enumerate(vars):
        meta.loc[meta['lobule'] == v,'domain'] = 1
    vars = ['VIII','IX','X'] #PZ/NZ
    for i,v in enumerate(vars):
        meta.loc[meta['lobule'] == v,'domain'] = 2
        
    return meta

def IPSC_reader():
    meta = IPSC_meta()
    os.chdir('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/mIPSCs/Combined/excel/clean') #Directory that contains IPSC !Excel! files
    output_df = pd.DataFrame()
    dirlist = os.listdir() 
    
    for filename in dirlist:
        if not filename.endswith('.xlsx'):
            continue 
        cellname = filename.split(".")[0]
        #print(meta.loc[meta.filename == filename.split(".")[0], "group"])
        print(cellname)
        zebrin_id = meta.loc[meta.filename == filename.split(".")[0],"group"].item()
        lobule = meta.loc[meta.filename == filename.split(".")[0],"lobule"].item()
        confirmed = meta.loc[meta.filename == filename.split(".")[0],"confirmed"].item()
        domain = meta.loc[meta.filename == filename.split(".")[0],"domain"].item()
        df = pd.read_excel(filename)
        df = df.drop(df.columns[0], axis=1)
        Frequency = len(df) / 120
        df['Zebrin']            = zebrin_id
        df['domain']            = domain
        df['Frequency (Hz)']    = Frequency
        df['Cellname']          = filename
        df['lobule']            = lobule
        df['confirmed']         = confirmed
        df['amplitude']         = df['amplitude'] * -1 #flips all negative values to positive values
        df['rise']              = df['rise'] * 1000 #converts seconds to miliseconds
        df['tau']               = df['tau'] * 1000 #converts seconds to miliseconds
        output_df = output_df.append(df, ignore_index = True)
    return(output_df)

def exclusion_list_ipsc(i):
    #creates a list of cells to exclude based on criteria
    exclusion_df = pd.DataFrame()
    
    #baseline
    baseline_cutoff = -500
    baseline_df     = i[i['baseline'] < baseline_cutoff]
    
    #amplitude
    amplitude_cutoff   = 120
    amplitude_df       = i[i['amplitude'] > amplitude_cutoff]
    
    #frequency
    frequency_cutoff1  = 15
    frequency_cutoff2  = 0.5
    frequency_df1      = i[i['Frequency (Hz)'] > frequency_cutoff1]
    frequency_df2      = i[i['Frequency (Hz)'] < frequency_cutoff2]
    
    #rise time
    rise_cutoff     = 8 
    rise_df         = i[i['rise'] > rise_cutoff]
    
    #tau
    tau_cutoff      = 20
    tau_df          = i[i['tau'] > tau_cutoff]
    
    #putting it all together
    exclusion_df = exclusion_df.append(baseline_df)
    exclusion_df = exclusion_df.append(amplitude_df)
    exclusion_df = exclusion_df.append(frequency_df1)
    exclusion_df = exclusion_df.append(frequency_df2)
    exclusion_df = exclusion_df.append(rise_df)
    exclusion_df = exclusion_df.append(tau_df)
    exclusion_df = exclusion_df.drop_duplicates()

    return exclusion_df

def export_IPSC_data():
    output_df = IPSC_reader()
    temp_median = output_df.groupby(['Zebrin','Cellname']).median().reset_index()
    temp_mean = output_df.groupby(['Zebrin','Cellname']).mean().reset_index()
    
    os.chdir('//storage.erasmusmc.nl/m/MyDocs/051071/My Documents/Desktop/Analysis_tools/output')
    temp_median.to_excel('./IPSC_mediandata.xlsx')
    temp_mean.to_excel('./IPSC_meandata.xlsx')
    return



############ PLOTTING ################
def plot_data_epsc():
    output_df    = EPSC_reader() 
    e = output_df.groupby(['Zebrin','Cellname']).median().reset_index()
    exclusion_df = exclusion_list_epsc(e)
    i_todrop     = exclusion_df.index
    e    = e.drop(i_todrop) 
    return e

def plot_data_ipsc():
    output_df = IPSC_reader()
    i = output_df.groupby(['Zebrin','Cellname']).median().reset_index()
    exclusion_df = exclusion_list_ipsc(i)
    i_todrop     = exclusion_df.index
    i    = i.drop(i_todrop) 
    return i

def plot_epsc_median(e):
    cm = 1/2.54
    #exclusion_df = exclusion_list_epsc(e)
    #i_todrop     = exclusion_df.index
    #plot_data    = e.drop(i_todrop) 
    plot_data    = e 
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    fig, ax = plt.subplots(1,4, figsize=(9.15*cm,3),dpi=300)

    for i, v in enumerate(vars):
    
        #sns.barplot(x = "Zebrin", y=v , data = plot_data, ax=ax[i], errwidth= 1.5,capsize= .25, errcolor = 'black',
        #            edgecolor='black', linewidth = 1, saturation = .5, order=['pos','neg'], ci='sd') 
        sns.boxplot(x='Zebrin',y=v, data=plot_data, ax=ax[i], order=['pos','neg'], fliersize=0, linewidth=0.5)
        sns.swarmplot(x="Zebrin", y=v, data= plot_data, ax=ax[i], dodge=True, palette = ['black','white'], edgecolor='black'
                      , linewidth=0.5, order=['pos','neg'], size=3, marker='^')        
        
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        
    #labels    
    ax[0].set_ylabel("Amplitude (pA)")
    ax[2].set_ylabel("Rise (ms)")
    ax[3].set_ylabel("Tau (ms)")
    ax[0].set_xlabel(" ")
    ax[1].set_xlabel(" ")
    ax[2].set_xlabel(" ")
    ax[3].set_xlabel(" ")
    ax[0].set_xticklabels(['+', '-'])
    ax[1].set_xticklabels(['+', '-'])
    ax[2].set_xticklabels(['+', '-'])
    ax[3].set_xticklabels(['+', '-'])

    #limits
    ax[0].set_ylim(0,50)
    ax[1].set_ylim(0,20)
    ax[2].set_ylim(0,10)
    ax[3].set_ylim(0,20)
    ax[0].set_yticks([0,10,20,30,40,50])
    ax[1].set_yticks([0,4,8,12,16,20])
    ax[2].set_yticks([0,2,4,6,8,10])
    ax[3].set_yticks([0,4,8,12,16,20])
    
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'mEPSC_boxplot.pdf'), dpi=300)
    return fig, ax

def plot_lobules(e):
    exclusion_df = exclusion_list_epsc(e)
    i_todrop     = exclusion_df.index
    plot_data    = e.drop(i_todrop)
    
    data = plot_data.loc[plot_data.confirmed == 'y']
    data = data.loc[data.lobule != '?']
    
    sns.set_style('ticks')
    sns.set_context('paper')
    fig, ax = plt.subplots(1,4)
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    for i,v in enumerate(vars):
        
        sns.barplot(x='domain',y=v,hue='Zebrin', data=data, ax=ax[i], edgecolor='black', capsize=.2
                    , hue_order=['pos','neg'])
        sns.stripplot(x='domain',y=v, hue='Zebrin',data=data, color='black', ax=ax[i], dodge=True
                    , hue_order=['pos','neg'])
    
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        
    return fig, ax

def plot_ipsc_median(i):
    cm=1/2.54
    #exclusion_df = exclusion_list_ipsc(i)
    #i_todrop     = exclusion_df.index
    #plot_data    = i.drop(i_todrop) 
    plot_data = i
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    fig, ax = plt.subplots(1,4,figsize= (9.15*cm,3),dpi=300)

    
    for i, v in enumerate(vars):
        #ax_x = int(i/2)
        #if i < 2:
        #    ax_y = i
        #else:
        #    ax_y = i-2
        sns.boxplot(x='Zebrin',y=v, data=plot_data, ax=ax[i], order=['pos','neg'], fliersize=0, linewidth=0.5)
        sns.swarmplot(x="Zebrin", y=v, data= plot_data, ax=ax[i], dodge=True, palette = ['black','white'], edgecolor='black'
                      , linewidth=0.5, order=['pos','neg'], size=3, marker='^')      
        
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        
    #labels    
    ax[0].set_ylabel("Amplitude (pA)")
    ax[2].set_ylabel("Rise (ms)")
    ax[3].set_ylabel("Tau (ms)")
    ax[0].set_xlabel(" ")
    ax[1].set_xlabel(" ")
    ax[2].set_xlabel(" ")
    ax[3].set_xlabel(" ")
    ax[0].set_xticklabels(['+', '-'])
    ax[1].set_xticklabels(['+', '-'])
    ax[2].set_xticklabels(['+', '-'])
    ax[3].set_xticklabels(['+', '-'])
    
    #limits
    ax[0].set_ylim(0,150)
    ax[1].set_ylim(0,20)
    ax[2].set_ylim(0,10)
    ax[3].set_ylim(0,20)
    ax[0].set_yticks([0,25,50,75,100,125,150])
    ax[1].set_yticks([0,4,8,12,16,20])
    ax[2].set_yticks([0,2,4,6,8,10])
    ax[3].set_yticks([0,4,8,12,16,20])
    
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc,'ipsc_boxplot.pdf'), dpi=300)
    return fig, ax

def complete_plot(e,i):
    cm=1/2.54
    exclusion_df = exclusion_list_ipsc(i)
    i_todrop     = exclusion_df.index
    plot_data_ipsc    = i.drop(i_todrop) 
    
    exclusion_df = exclusion_list_epsc(e)
    i_todrop     = exclusion_df.index
    plot_data_epsc    = e.drop(i_todrop) 
    
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    fig, ax = plt.subplots(2,4,figsize= (9.15*cm,6),dpi=300)
    
    for i, v in enumerate(vars):
        sns.boxplot(x='Zebrin',y=v, data=plot_data_epsc, ax=ax[0,i], order=['pos','neg'], fliersize=0, linewidth=0.5)
        sns.stripplot(x="Zebrin", y=v, data= plot_data_epsc, ax=ax[0,i], dodge=True, palette = ['black','white'], edgecolor='black'
                      ,linewidth=0.5, order=['pos','neg'], size=3, marker='^')      
        sns.boxplot(x='Zebrin',y=v, data=plot_data_ipsc, ax=ax[1,i], order=['pos','neg'], fliersize=0, linewidth=0.5)
        sns.stripplot(x="Zebrin", y=v, data= plot_data_ipsc, ax=ax[1,i], dodge=True, palette = ['black','white'], edgecolor='black'
                      ,linewidth=0.5, order=['pos','neg'], size=3, marker='^')     
        
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.set_xticklabels(['+','-'])

        
    ax[0,0].set_ylim(0,50)
    ax[0,1].set_ylim(0,20)
    ax[0,2].set_ylim(0,10)
    ax[0,3].set_ylim(0,20)
    ax[1,0].set_ylim(0,150)
    ax[1,1].set_ylim(0,20)
    ax[1,2].set_ylim(0,10)
    ax[1,3].set_ylim(0,20)
    ax[1,0].set_yticks([0,25,50,75,100,125,150])
    ax[0,1].set_yticks([0,4,8,12,16,20])
    ax[1,1].set_yticks([0,4,8,12,16,20])
    ax[0,3].set_yticks([0,4,8,12,16,20])
    ax[1,3].set_yticks([0,4,8,12,16,20])
    ax[0,0].set_xlabel(' ')
    ax[0,1].set_xlabel(' ')
    ax[0,2].set_xlabel(' ')
    ax[0,3].set_xlabel(' ')
    ax[1,0].set_xlabel('EAAT4')
    ax[1,1].set_xlabel('EAAT4')
    ax[1,2].set_xlabel('EAAT4')
    ax[1,3].set_xlabel('EAAT4')
    ax[0,0].set_ylabel("Amplitude (pA)")
    ax[0,2].set_ylabel("Rise (ms)")
    ax[0,3].set_ylabel("Tau (ms)")
    ax[1,0].set_ylabel("Amplitude (pA)")
    ax[1,2].set_ylabel("Rise (ms)")
    ax[1,3].set_ylabel("Tau (ms)")
    
    fig.tight_layout()
    fig.align_ylabels()
    fig.savefig(os.path.join(output_loc,'complete_psc_plot.pdf'), dpi=300)

################################################## STATISTICS ####################################
def medstats_epsc(e):
    exclusion_df = exclusion_list_epsc(e)
    i_todrop     = exclusion_df.index
    e            = e.drop(i_todrop)
    
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    for v in vars:
        print(v, tstd(e.loc[e.Zebrin == "pos", v]), 'Z+ standard deviation')
        print(v, tstd(e.loc[e.Zebrin == "neg", v]), 'Z- standard deviation')
        print(v, levene(e.loc[e.Zebrin == "pos", v], e.loc[e.Zebrin == "neg", v]))
        print(v, ttest_ind(e.loc[e.Zebrin == "pos", v], e.loc[e.Zebrin == "neg", v]), "\n")
        
def medstats_ipsc(i):
    exclusion_df = exclusion_list_ipsc(i)
    i_todrop     = exclusion_df.index
    i            = i.drop(i_todrop)
    
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    for v in vars:
        print(v, tstd(i.loc[i.Zebrin == "pos", v]), 'Z+ standard deviation')
        print(v, tstd(i.loc[i.Zebrin == "neg", v]), 'Z- standard deviation')
        print(v, levene(i.loc[i.Zebrin == "pos", v], i.loc[i.Zebrin == "neg", v]))
        print(v, ttest_ind(i.loc[i.Zebrin == "pos", v], i.loc[i.Zebrin == "neg", v]), "\n")   
        
        
        
### Just messing around stuff 
def epsc_violin_plot(e):
    fig, ax = plt.subplots(1,4)
    vars = ['amplitude','Frequency (Hz)','rise','tau']
    for i,v in enumerate(vars):
        sns.violinplot(data=e,x='Zebrin',y=v,ax=ax[i])
    fig.tight_layout()
    
    return fig, ax
                
def ipsc_violin_plot(i):
    fig, ax = plt.subplots(1,4)
    vars = ['amplitude','Frequency (Hz)','rise','tau']
    for c,v in enumerate(vars):
        sns.violinplot(data=i,x='Zebrin',y=v,ax=ax[c])
    fig.tight_layout()
    
    return fig, ax 

def lobule_plot(e,i):
    cm = 1/2.54
    #epsc_data = e.loc[e.confirmed == 'y']
    #ipsc_data = i.loc[i.confirmed == 'y']
    e.dropna()
    i.dropna()
    epsc_data = e
    ipsc_data = i
    lobules = [0,1,2]
    
    fig, ax = plt.subplots(2,4, figsize = ((18.3*cm), 3))
    vars = ['amplitude','Frequency (Hz)','rise','tau']
    for i,v in enumerate(vars):
        sns.boxplot(x='domain',y=v,hue='Zebrin', data=epsc_data, ax=ax[0,i], order=lobules, hue_order=['pos','neg'], linewidth=0.5)
        sns.boxplot(x='domain',y=v,hue='Zebrin', data=ipsc_data, ax=ax[1,i], order=lobules, hue_order=['pos','neg'], linewidth=0.5)
        sns.stripplot(x='domain',y=v,hue='Zebrin', data=epsc_data, ax=ax[0,i], order=lobules, hue_order=['pos','neg'], dodge=True,
                      color = 'white', edgecolor='black', linewidth=0.5, jitter=0.05, marker='^')
        sns.stripplot(x='domain',y=v,hue='Zebrin', data=ipsc_data, ax=ax[1,i], order=lobules, hue_order=['pos','neg'], dodge=True,
                      color = 'white', edgecolor='black', linewidth=0.5, jitter=0.05, marker='^')
        
        for a in ax.flatten():
            a.spines['top'].set_visible(False)
            a.spines['right'].set_visible(False)
            a.legend().set_visible(False)
            a.set_xticklabels(['AZ','CZ','PZ/NZ'])

        ax[0,0].set_xlabel(' ')
        ax[0,1].set_xlabel(' ')
        ax[0,2].set_xlabel(' ')
        ax[0,3].set_xlabel(' ')
        ax[0,0].set_ylabel('Amplitude (pA)')
        ax[0,2].set_ylabel('Rise time (ms)')
        ax[0,3].set_ylabel('Tau (ms)')
        ax[1,0].set_xlabel('Zone')
        ax[1,1].set_xlabel('Zone')
        ax[1,2].set_xlabel('Zone')
        ax[1,3].set_xlabel('Zone')
        ax[1,0].set_ylabel('Amplitude (pA)')
        ax[1,2].set_ylabel('Rise time (ms)')
        ax[1,3].set_ylabel('Tau (ms)')
        
        ax[0,0].set_ylim(0,50)
        ax[0,0].set_yticks([0,10,20,30,40,50])
        ax[0,1].set_ylim(0,10)
        ax[0,2].set_ylim(0,10)
        ax[0,3].set_ylim(0,20)
        ax[1,0].set_ylim(0,250)
        ax[1,0].set_yticks([0,50,100,150,200,250])
        ax[1,1].set_ylim(0,10)
        ax[1,2].set_ylim(0,10)
        ax[1,3].set_ylim(0,20)
    fig.tight_layout()
    fig.align_ylabels()
    fig.savefig(os.path.join(output_loc, "minidomains.pdf"))
    return fig, ax

#needs all event data, so not meaned
def distribution(epsc_data, ipsc_data):
    cm = 1/2.54
    vars = ['amplitude','rise','tau']
    
    for v in vars:
        epsc_data_norm =  (epsc_data[v] - epsc_data[v].min()) / (epsc_data[v].max() - epsc_data[v].min()) 
        ipsc_data_norm =  (ipsc_data[v] - ipsc_data[v].min()) / (ipsc_data[v].max() - ipsc_data[v].min())
        
    fig, ax = plt.subplots(2,3, figsize=(9.15*cm,3))

    for i,v in enumerate(vars):
        sns.kdeplot(data=ipsc_data,x=v,ax=ax[1,i], hue='Zebrin', hue_order=['pos','neg'], 
                    fill=False, multiple='layer')
        sns.kdeplot(data=epsc_data,x=v,ax=ax[0,i], hue='Zebrin', hue_order=['pos','neg'],  
                    fill=False, multiple='layer')
        
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.get_legend().remove()
        a.set_ylabel(' ')
        a.set_yticklabels([' '])
        a.tick_params(axis='y',which='both',bottom=True,top=False,labelbottom=True, left=False)
    
    for i in range(3):
        ax[0,i].set_xlabel(' ')
    ax[0,0].set_xlabel(' ')
    ax[0,1].set_xlabel(' ')
    ax[0,2].set_xlabel(' ')
    ax[1,0].set_xlabel('Amplitude (pA)')
    ax[1,1].set_xlabel('Rise (ms)')
    ax[1,2].set_xlabel('Tau (ms)')
    ax[0,0].set_ylabel('Density')
    ax[1,0].set_ylabel('Density')
    # ax[0,0].set_ylim(0,0.06)
    # ax[0,1].set_ylim(0,0.15)
    # ax[0,2].set_ylim(0,0.05)
    # ax[1,0].set_ylim(0,0.008)
    # ax[1,1].set_ylim(0,0.22)
    # ax[1,2].set_ylim(0,0.05)
    ax[0,0].set_xlim(0,100)
    ax[0,1].set_xlim(0,20)
    ax[0,2].set_xlim(0,50)
    ax[1,0].set_xlim(0,250)
    ax[1,1].set_xlim(0,20)
    ax[1,2].set_xlim(0,100)

    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "minidistributions.pdf"))
    return fig, ax

def compare_confirmed(e, i): 
    #Mainly to check if the distributions of the confirmed populations are similar to the whole population
    #Or bluntly, if I am any good at actually determining EAAT4 ID in vitro
    epsc = e.loc[e.confirmed != 'n']
    ipsc = i.loc[i.confirmed != 'n']
    
    vars = ['amplitude','Frequency (Hz)','rise','tau']
    fig, ax = plt.subplots(4,4, figsize=(6,6))
    for x, v in enumerate(vars):
        sns.boxplot(x='Zebrin',y=v,data=epsc,ax=ax[0,x])
        sns.boxplot(x='Zebrin',y=v,data=e,ax=ax[1,x])
    
        sns.boxplot(x='Zebrin',y=v,data=ipsc,ax=ax[2,x])
        sns.boxplot(x='Zebrin',y=v,data=i,ax=ax[3,x])
        
    fig.tight_layout()
    return fig, ax

def bar_plot_epsc(e):
    cm=1/2.54
    plot_data = e
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    fig, ax = plt.subplots(1,4,figsize= (9.15*cm,3),dpi=300)
    for i, v in enumerate(vars):
        sns.barplot(x='Zebrin',y=v, data=plot_data, ax=ax[i], order=['pos','neg'], errorbar=('ci',95), width=0.5, errcolor='black', errwidth=0.5, 
                    capsize=0.3, facecolor=(0,0,0,0), edgecolor=(["#5BE12D","#EC008C"]), lw=0.5)
        sns.stripplot(x="Zebrin", y=v, data= plot_data, ax=ax[i], dodge=False, palette = ['black','white'], edgecolor='black'
                      , linewidth=0.5,order=['pos','neg'], size=2, marker='^', jitter=0.05)      
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
    #labels    
    ax[0].set_ylabel("Amplitude (pA)")
    ax[2].set_ylabel("Rise (ms)")
    ax[3].set_ylabel("Tau (ms)")
    ax[0].set_xlabel(" ")
    ax[1].set_xlabel(" ")
    ax[2].set_xlabel(" ")
    ax[3].set_xlabel(" ")
    ax[0].set_xticklabels(['+', '-'])
    ax[1].set_xticklabels(['+', '-'])
    ax[2].set_xticklabels(['+', '-'])
    ax[3].set_xticklabels(['+', '-'])

    #limits
    ax[0].set_ylim(0,50)
    ax[1].set_ylim(0,20)
    ax[2].set_ylim(0,10)
    ax[3].set_ylim(0,20)
    ax[0].set_yticks([0,10,20,30,40,50])
    ax[1].set_yticks([0,4,8,12,16,20])
    ax[2].set_yticks([0,2,4,6,8,10])
    ax[3].set_yticks([0,4,8,12,16,20])
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc,'epsc_boxplot.png'), dpi=1000)
    return fig, ax

def bar_plot_ipsc(i):
    cm=1/2.54
    plot_data = i
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    fig, ax = plt.subplots(1,4,figsize= (9.15*cm,3),dpi=300)
    for i, v in enumerate(vars):
        sns.barplot(x='Zebrin',y=v, data=plot_data, ax=ax[i], order=['pos','neg'], errorbar=('ci',95), width=0.5, errcolor='black', errwidth=0.5, 
                    capsize=0.3, facecolor=(0,0,0,0), edgecolor=(["#5BE12D","#EC008C"]), lw=0.5)
        sns.stripplot(x="Zebrin", y=v, data= plot_data, ax=ax[i], dodge=False, palette = ['black','white'], edgecolor='black'
                      , linewidth=0.5,order=['pos','neg'], size=2, marker='^', jitter=0.05)      
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
    #labels    
    ax[0].set_ylabel("Amplitude (pA)")
    ax[2].set_ylabel("Rise (ms)")
    ax[3].set_ylabel("Tau (ms)")
    ax[0].set_xlabel(" ")
    ax[1].set_xlabel(" ")
    ax[2].set_xlabel(" ")
    ax[3].set_xlabel(" ")
    ax[0].set_xticklabels(['+', '-'])
    ax[1].set_xticklabels(['+', '-'])
    ax[2].set_xticklabels(['+', '-'])
    ax[3].set_xticklabels(['+', '-'])
    
    #limits
    ax[0].set_ylim(0,150)
    ax[1].set_ylim(0,20)
    ax[2].set_ylim(0,10)
    ax[3].set_ylim(0,20)
    ax[0].set_yticks([0,25,50,75,100,125,150])
    ax[1].set_yticks([0,4,8,12,16,20])
    ax[2].set_yticks([0,2,4,6,8,10])
    ax[3].set_yticks([0,4,8,12,16,20])
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc,'ipsc_boxplot.png'), dpi=1000)
    return fig, ax

def extract(data):
    data.to_excel(os.path.join(output_loc,'excel.xlsx'))
    return
    
if __name__ == "__main__":
    base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
    output_loc = os.path.join(base_folder, "Analysis_tools/output/")
    plt.style.use('seaborn-paper')
    color = ["#5BE12D","#EC008C"] #Original green/red
    sns.set_palette(sns.color_palette(color))
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['xtick.minor.width'] = 0.3
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['ytick.minor.width'] = 0.3
    fonts = {"font.size":7, "axes.labelsize":7, "ytick.labelsize":7, "xtick.labelsize":7}
    plt.rcParams.update(fonts)
