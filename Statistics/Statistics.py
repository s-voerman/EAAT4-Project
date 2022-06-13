import matplotlib as mpl
import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
from scipy.stats import levene, ttest_ind, tstd, sem

def descriptives(data, output_loc):
    meta = {}
    data = data.loc[data.time > -6]
    plot_data  = data.groupby(['group','cell','pre_post']).median().reset_index()
    temp_data  = plot_data.groupby('pre_post')
    pre_data   = temp_data.get_group('pre') 
    post_data  = temp_data.get_group('post')
    pre_pos = pre_data.loc[pre_data.group == 'pos'] 
    pre_neg = pre_data.loc[pre_data.group == 'neg'] 
    post_pos = post_data.loc[post_data.group == 'pos']
    post_neg = post_data.loc[post_data.group == 'neg']
    
    str_list = ['pre_pos', 'post_pos', 'pre_neg', 'post_neg'] 
    list = [pre_pos, post_pos, pre_neg, post_neg] 

    vars = ['norm_EPSC1_amplitude'] #Still would like to add multiple variable support!
    meta = {}
    
    for c,i in enumerate(list):
        for v in vars:
            mean = np.mean(i.loc[:,v])
            sderr = sem(i.loc[:,v])
            sd = tstd(i.loc[:,v])
    
            meta[str_list[c],v,'Mean'] =  mean
            meta[str_list[c],v,'SEM'] = sderr
            meta[str_list[c],v,'SD'] = sd
            
    stat_df = pd.DataFrame([meta])
    stat_df.to_excel(os.path.join(output_loc, 'LTP_PT ' + str(v) + '.xlsx'))
    return stat_df

def plasticity_stats(data, output_loc):
    data = data.loc[data.time > -6]
    plot_data  = data.groupby(['group','cell','pre_post']).median().reset_index()
    temp_data  = plot_data.groupby('pre_post')
    pre_data   = temp_data.get_group('pre')
    post_data  = temp_data.get_group('post')
    
    vars = ["norm_EPSC1_amplitude"]
    for v in vars:
        print(v, levene(pre_data.loc[pre_data.group == "pos", v], post_data.loc[post_data.group == "pos", v]), 'prepost_pos/pos')
        print(v, ttest_ind(pre_data.loc[pre_data.group == "pos", v], post_data.loc[post_data.group == "pos", v]), "\n")
    for v in vars:
        print(v, levene(pre_data.loc[pre_data.group == "neg", v], post_data.loc[post_data.group == "neg", v]), 'prepost_neg/neg')
        print(v, ttest_ind(pre_data.loc[pre_data.group == "neg", v], post_data.loc[post_data.group == "neg", v]), "\n")
    for v in vars:
        print(v, levene(pre_data.loc[pre_data.group == "pos", v], pre_data.loc[pre_data.group == "neg", v]), 'pre_neg/pos')
        print(v, ttest_ind(pre_data.loc[pre_data.group == "pos", v], pre_data.loc[pre_data.group == "neg", v]), "\n")
    for v in vars:
        print(v, levene(post_data.loc[post_data.group == "pos", v], post_data.loc[post_data.group == "neg", v]), 'post_neg/pos')
        print(v, ttest_ind(post_data.loc[post_data.group == "pos", v], post_data.loc[post_data.group == "neg", v]), "\n")
    return

def descriptives_25(data, output_loc):
    meta = {}
    data = data.loc[(data.time > 19) & (data.time < 26)]  #used for stats @ 20-25 mins, comment out all /// for this
    plot_data  = data.groupby(['group','cell','pre_post']).median().reset_index()
    temp_data  = plot_data.groupby('pre_post')
    post_data  = temp_data.get_group('post')
    post_pos = post_data.loc[post_data.group == 'pos']
    post_neg = post_data.loc[post_data.group == 'neg']

    str_list = ['post_pos', 'post_neg']
    list = [post_pos, post_neg]
    vars = ['norm_EPSC1_amplitude'] #Still would like to add multiple variable support!
    meta = {}
    
    for c,i in enumerate(list):
        for v in vars:
            mean = np.mean(i.loc[:,v])
            sderr = sem(i.loc[:,v])
            sd = tstd(i.loc[:,v])
    
            meta[str_list[c],v,'Mean'] =  mean
            meta[str_list[c],v,'SEM'] = sderr
            meta[str_list[c],v,'SD'] = sd
            
    stat_df = pd.DataFrame([meta])
    stat_df.to_excel(os.path.join(output_loc, 'LTD ' + str(v) + ' ' + str(25) + '.xlsx'))
    return stat_df

def plasticity_stats_25(data, output_loc):
    data = data.loc[(data.time > 19) & (data.time < 26)]
    plot_data  = data.groupby(['group','cell','pre_post']).median().reset_index()
    temp_data  = plot_data.groupby('pre_post')
    post_data  = temp_data.get_group('post')
    
    vars = ["norm_EPSC1_amplitude"]
    for v in vars:
        print(v, levene(post_data.loc[post_data.group == "pos", v], post_data.loc[post_data.group == "neg", v]), 'post_neg/pos')
        print(v, ttest_ind(post_data.loc[post_data.group == "pos", v], post_data.loc[post_data.group == "neg", v]), "\n")
    return

def descriptives_zones(data, output_loc):
    meta = {}
    data = data.loc[(data.time > 19) & (data.time < 26)]  #used for stats @ 20-25 mins, comment out all /// for this
    plot_data  = data.groupby(['group','cell','pre_post','domain']).median().reset_index()
    
    AZ_data_pos = plot_data.loc[(plot_data.domain == 0) & (plot_data.group == 'pos')]
    AZ_data_neg = plot_data.loc[(plot_data.domain == 0) & (plot_data.group == 'neg')]
    CZ_data_pos = plot_data.loc[(plot_data.domain == 1) & (plot_data.group == 'pos')]
    CZ_data_neg = plot_data.loc[(plot_data.domain == 1) & (plot_data.group == 'neg')]
    PZ_data_pos = plot_data.loc[(plot_data.domain == 2) & (plot_data.group == 'pos')]
    PZ_data_neg = plot_data.loc[(plot_data.domain == 2) & (plot_data.group == 'neg')]
    
    str_list = ['AZ_data_pos', 'AZ_data_neg', 'CZ_data_pos', 'CZ_data_neg', 'PZ_data_pos','PZ_data_neg']
    list = [AZ_data_pos, AZ_data_neg, CZ_data_pos, CZ_data_neg, PZ_data_pos, PZ_data_neg]
    vars = ['norm_EPSC1_amplitude'] #Still would like to add multiple variable support!
    meta = {}
    
    for c,i in enumerate(list):
        for v in vars:
            mean = np.mean(i.loc[:,v])
            sderr = sem(i.loc[:,v])
            sd = tstd(i.loc[:,v])
    
            meta[str_list[c],v,'Mean'] =  mean
            meta[str_list[c],v,'SEM'] = sderr
            meta[str_list[c],v,'SD'] = sd
            
    stat_df = pd.DataFrame([meta])
    stat_df.to_excel(os.path.join(output_loc, 'LTP_RT ' + str(v) + ' ' + 'zones' + '.xlsx'))
    return stat_df

def plasticity_stats_zones(data, output_loc):
    data = data.loc[(data.time > 19) & (data.time < 26)]
    data = data.loc[(data.time > 19) & (data.time < 26)]  #used for stats @ 20-25 mins, comment out all /// for this
    plot_data  = data.groupby(['group','cell','pre_post','domain']).median().reset_index()
    
    AZ_data = plot_data.loc[plot_data.domain == 0]
    CZ_data = plot_data.loc[plot_data.domain == 1]
    PZ_data = plot_data.loc[plot_data.domain == 2]
    
    str_list = ['AZ_data_pos', 'AZ_data_neg', 'CZ_data_pos', 'CZ_data_neg', 'PZ_data_pos','PZ_data_neg']
    list = [AZ_data_pos, AZ_data_neg, CZ_data_pos, CZ_data_neg, PZ_data_pos, PZ_data_neg]
    
    vars = ["norm_EPSC1_amplitude"]
    for v in vars:
        print(v, levene(AZ_data.loc[AZ_data.group == "pos", v], AZ_data.loc[AZ_data.group == "neg", v]), 'post_neg/pos AZ')
        print(v, ttest_ind(AZ_data.loc[AZ_data.group == "pos", v], AZ_data.loc[AZ_data.group == "neg", v]), "\n")
    for v in vars:
        print(v, levene(CZ_data.loc[CZ_data.group == "pos", v], CZ_data.loc[CZ_data.group == "neg", v]), 'post_neg/pos CZ')
        print(v, ttest_ind(CZ_data.loc[CZ_data.group == "pos", v], CZ_data.loc[CZ_data.group == "neg", v]), "\n")
    for v in vars:
        print(v, levene(PZ_data.loc[PZ_data.group == "pos", v], PZ_data.loc[PZ_data.group == "neg", v]), 'post_neg/pos CZ')
        print(v, ttest_ind(PZ_data.loc[PZ_data.group == "pos", v], PZ_data.loc[PZ_data.group == "neg", v]), "\n")

    return

def descriptives_epscs(e, output_loc):
    meta = {}
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    epsc_pos = e.loc[e.Zebrin == 'pos']
    epsc_neg = e.loc[e.Zebrin == 'neg']
    str_list = ['epsc_pos','epsc_neg']
    list = [epsc_pos, epsc_neg]
    for c,j in enumerate(list):
        for v in vars:
            mean = np.mean(j.loc[:,v])
            sderr = sem(j.loc[:,v])
            sd = tstd(j.loc[:,v])
    
            meta[str_list[c],v,'Mean'] =  mean
            meta[str_list[c],v,'SEM'] = sderr
            meta[str_list[c],v,'SD'] = sd
    stat_df = pd.DataFrame([meta])
    stat_df.to_excel(os.path.join(output_loc, 'epsc ' + str(v) + '.xlsx'))
    return

def epsc_stats_zones(e, output_loc):
    AZ_data = e.loc[e.domain == 0]
    CZ_data = e.loc[e.domain == 1]
    PZ_data = e.loc[e.domain == 2]
        
    vars = ["Frequency (Hz)","amplitude", "rise", "tau"]
    for v in vars:
        print(v, levene(AZ_data.loc[AZ_data.Zebrin == "pos", v], AZ_data.loc[AZ_data.Zebrin == "neg", v]), 'AZ')
        print(v, ttest_ind(AZ_data.loc[AZ_data.Zebrin == "pos", v], AZ_data.loc[AZ_data.Zebrin == "neg", v]), "\n")
    for v in vars:
        print(v, levene(CZ_data.loc[CZ_data.Zebrin == "pos", v], CZ_data.loc[CZ_data.Zebrin == "neg", v]), 'CZ')
        print(v, ttest_ind(CZ_data.loc[CZ_data.Zebrin == "pos", v], CZ_data.loc[CZ_data.Zebrin == "neg", v]), "\n")
    for v in vars:
        print(v, levene(PZ_data.loc[PZ_data.Zebrin == "pos", v], PZ_data.loc[PZ_data.Zebrin == "neg", v]), 'CZ')
        print(v, ttest_ind(PZ_data.loc[PZ_data.Zebrin == "pos", v], PZ_data.loc[PZ_data.Zebrin == "neg", v]), "\n")

    return

def descriptives_epscs_zones(e, output_loc):
    meta = {}
    AZ_data_pos = e.loc[(e.domain == 0) & (e.Zebrin == 'pos')]
    AZ_data_neg = e.loc[(e.domain == 0) & (e.Zebrin == 'neg')]
    CZ_data_pos = e.loc[(e.domain == 1) & (e.Zebrin == 'pos')]
    CZ_data_neg = e.loc[(e.domain == 1) & (e.Zebrin == 'neg')]
    PZ_data_pos = e.loc[(e.domain == 2) & (e.Zebrin == 'pos')]
    PZ_data_neg = e.loc[(e.domain == 2) & (e.Zebrin == 'neg')]
    
    str_list = ['AZ_data_pos', 'AZ_data_neg', 'CZ_data_pos', 'CZ_data_neg', 'PZ_data_pos','PZ_data_neg']
    list = [AZ_data_pos, AZ_data_neg, CZ_data_pos, CZ_data_neg, PZ_data_pos, PZ_data_neg]
    
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    for c,i in enumerate(list):
        for v in vars:
            mean = np.mean(i.loc[:,v])
            sderr = sem(i.loc[:,v])
            sd = tstd(i.loc[:,v])
    
            meta[str_list[c],v,'Mean'] =  mean
            meta[str_list[c],v,'SEM'] = sderr
            meta[str_list[c],v,'SD'] = sd
            
    stat_df = pd.DataFrame([meta])
    stat_df.to_excel(os.path.join(output_loc, 'epsc ' + str(v) + ' ' + 'zones' + '.xlsx'))
    return stat_df

def descriptives_ipscs(i, output_loc):
    meta = {}
    vars = ["Frequency (Hz)","amplitude", "rise", "tau"]
    ipsc_pos = i.loc[i.Zebrin == 'pos']
    ipsc_neg = i.loc[i.Zebrin == 'neg']
    str_list = ['ipsc_pos','ipsc_neg']
    list = [ipsc_pos, ipsc_neg]
    for c,j in enumerate(list):
        for v in vars:
            mean = np.mean(j.loc[:,v])
            sderr = sem(j.loc[:,v])
            sd = tstd(j.loc[:,v])
    
            meta[str_list[c],v,'Mean'] =  mean
            meta[str_list[c],v,'SEM'] = sderr
            meta[str_list[c],v,'SD'] = sd
    stat_df = pd.DataFrame([meta])
    stat_df.to_excel(os.path.join(output_loc, 'ipsc ' + str(v) + '.xlsx'))
    return

def ipsc_stats_zones(i, output_loc):
    AZ_data = i.loc[i.domain == 0]
    CZ_data = i.loc[i.domain == 1]
    PZ_data = i.loc[i.domain == 2]
        
    vars = ["Frequency (Hz)","amplitude", "rise", "tau"]
    for v in vars:
        print(v, levene(AZ_data.loc[AZ_data.Zebrin == "pos", v], AZ_data.loc[AZ_data.Zebrin == "neg", v]), 'AZ')
        print(v, ttest_ind(AZ_data.loc[AZ_data.Zebrin == "pos", v], AZ_data.loc[AZ_data.Zebrin == "neg", v]), "\n")
    for v in vars:
        print(v, levene(CZ_data.loc[CZ_data.Zebrin == "pos", v], CZ_data.loc[CZ_data.Zebrin == "neg", v]), 'CZ')
        print(v, ttest_ind(CZ_data.loc[CZ_data.Zebrin == "pos", v], CZ_data.loc[CZ_data.Zebrin == "neg", v]), "\n")
    for v in vars:
        print(v, levene(PZ_data.loc[PZ_data.Zebrin == "pos", v], PZ_data.loc[PZ_data.Zebrin == "neg", v]), 'CZ')
        print(v, ttest_ind(PZ_data.loc[PZ_data.Zebrin == "pos", v], PZ_data.loc[PZ_data.Zebrin == "neg", v]), "\n")

    return

def descriptives_ipscs_zones(i, output_loc):
    meta = {}
    AZ_data_pos = i.loc[(i.domain == 0) & (i.Zebrin == 'pos')]
    AZ_data_neg = i.loc[(i.domain == 0) & (i.Zebrin == 'neg')]
    CZ_data_pos = i.loc[(i.domain == 1) & (i.Zebrin == 'pos')]
    CZ_data_neg = i.loc[(i.domain == 1) & (i.Zebrin == 'neg')]
    PZ_data_pos = i.loc[(i.domain == 2) & (i.Zebrin == 'pos')]
    PZ_data_neg = i.loc[(i.domain == 2) & (i.Zebrin == 'neg')]
    
    str_list = ['AZ_data_pos', 'AZ_data_neg', 'CZ_data_pos', 'CZ_data_neg', 'PZ_data_pos','PZ_data_neg']
    list = [AZ_data_pos, AZ_data_neg, CZ_data_pos, CZ_data_neg, PZ_data_pos, PZ_data_neg]
    
    vars = ["Frequency (Hz)","amplitude", "rise", "tau"]
    for c,j in enumerate(list):
        for v in vars:
            mean = np.mean(j.loc[:,v])
            sderr = sem(j.loc[:,v])
            sd = tstd(j.loc[:,v])
    
            meta[str_list[c],v,'Mean'] =  mean
            meta[str_list[c],v,'SEM'] = sderr
            meta[str_list[c],v,'SD'] = sd
            
    stat_df = pd.DataFrame([meta])
    stat_df.to_excel(os.path.join(output_loc, 'ipsc ' + str(v) + ' ' + 'zones' + '.xlsx'))
    return stat_df

def mini_seal_descriptives(output_loc):
    meta = {}
    seal_data = pd.read_excel('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/seal_data_combined.xlsx')
    pos_data = seal_data.loc[seal_data.group == 'pos']
    neg_data = seal_data.loc[seal_data.group == 'neg']
    
    str_list = ['pos_data','neg_data']
    list = [pos_data, neg_data]
    vars = ["Rm","Cm","Tau"]
    for c,j in enumerate(list):
        for v in vars:
            mean = np.mean(j.loc[:,v])
            sderr = sem(j.loc[:,v])
            sd = tstd(j.loc[:,v])
    
            meta[str_list[c],v,'Mean'] =  mean
            meta[str_list[c],v,'SEM'] = sderr
            meta[str_list[c],v,'SD'] = sd
            
    stat_df = pd.DataFrame([meta])
    stat_df.to_excel(os.path.join(output_loc, 'ipsc ' + str(v) + ' ' + 'zones' + '.xlsx'))
    return stat_df

def medstats_epsc(e):
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    for v in vars:
        print(v, levene(e.loc[e.Zebrin == "pos", v], e.loc[e.Zebrin == "neg", v]))
        print(v, ttest_ind(e.loc[e.Zebrin == "pos", v], e.loc[e.Zebrin == "neg", v]), "\n")
        
def medstats_ipsc(i):
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    for v in vars:
        print(v, levene(i.loc[i.Zebrin == "pos", v], i.loc[i.Zebrin == "neg", v]))
        print(v, ttest_ind(i.loc[i.Zebrin == "pos", v], i.loc[i.Zebrin == "neg", v]), "\n") 

if __name__ == "__main__":
    base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
    data_loc = os.path.join(base_folder, "Data/LTD/Data")
    meta_loc = os.path.join(base_folder, "Metadata/ltdmeta_select.xlsx")
    output_loc = os.path.join(base_folder, "Analysis_tools/output/")
