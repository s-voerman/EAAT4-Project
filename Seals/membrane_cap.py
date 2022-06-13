import numpy as np
import pandas as pd
import os
from neo.rawio import AxonRawIO
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.stats import levene, ttest_ind
from scipy.stats import tstd

def read_abf(file_loc):
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header["signal_channels"][0][2])
    data = a.get_analogsignal_chunk() * a.header["signal_channels"][0][5]
    return data[:, 0], sr

def read_masterfile():
    meta_loc = '//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Metadata'
    os.chdir(meta_loc)
    ipsc_df = pd.read_excel('mipsc_metafile.xlsx') 
    epsc_df = pd.read_excel('mepsc_metafile.xlsx')
    for i,x in ipsc_df.iterrows():
        file_loc = x['name'] + ' ' + str(x.cell) + '-' + str(x.seal)
        file_loc = file_loc + '.ABF'
        ipsc_df.loc[i,'file_loc'] = file_loc
    for i,x in epsc_df.iterrows():
        file_loc = x['name'] + ' ' + str(x.cell) + '-' + str(x.seal)
        file_loc = file_loc + '.ABF'
        epsc_df.loc[i,'file_loc'] = file_loc
    return ipsc_df, epsc_df

def test():
    ipsc_df, epsc_df = read_masterfile()
    
    ipsc_loc = '//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/mIPSCs/Combined/seals/abf'
    epsc_loc = '//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/mEPSCs/Combined/seals/abf'
    
    sr=50000
    baseline_start   = int(0.65 * sr)
    baseline_end     = int(0.7 * sr)
    peak_start       = int(0.7 * sr)
    peak_end         = int(0.725 * sr)
    inp_start        = int(0.780 * sr)
    inp_end          = int(0.8 * sr)
    cm_peak_start    = int(0.8*sr)
    cm_base_start    = int(0.9 * sr)
    
    meta = {}
    output_df = pd.DataFrame()
    
    os.chdir(ipsc_loc)
    for i,x in ipsc_df.iterrows():
        Ri, Rs, Rm, tau, Cm = [], [], [],  [], []
        data, sr = read_abf(x.file_loc)
        for start in np.arange(0, len(data), sr):
            baseline = np.mean(data[start + baseline_start:start + baseline_end])
            peak     = np.min(data[start + peak_start:start + peak_end])
            inp      = np.mean(data[start + inp_start:start + inp_end])
            
            input_resistance = -0.01/((inp - baseline)*1e-12)
            series_resistance = -0.01/((peak - baseline)*1e-12)
            membrane_resistance = input_resistance - series_resistance 
            input_resistance = input_resistance / 1e6
            membrane_resistance = membrane_resistance / 1e6
            series_resistance = series_resistance / 1e6
            
            tau_perc = 1 - (1-1/np.exp(1))
            y_data = data[start+cm_peak_start:start + cm_base_start]
            x_data = np.arange(0,len(y_data),1)
            temp = np.max(y_data) - np.mean(data[start + int(0.9*sr):start+int(0.95*sr)])
            Itau = tau_perc * temp
            y_data2 = np.full(len(x_data),Itau)
            
            plt.plot(x_data,y_data,'-')
            plt.plot(x_data,y_data2,'-')
            idx = np.argwhere(np.diff(np.sign(y_data - y_data2))).flatten()
            plt.plot(x_data[idx],y_data[idx],'ro')
            try:
                Tau = idx[len(idx)-1] * 1000/sr        
            except:
                pass
            cm = (Tau / 1000)/ (membrane_resistance * 1e6) * 1e12
            
            Ri.append(input_resistance)
            Rs.append(series_resistance)
            Rm.append(membrane_resistance)
            tau.append(Tau)
            Cm.append(cm)
            meta["Ri"] = np.mean(Ri)
            meta["Rs"] = np.mean(Rs)
            meta["Rm"] = np.mean(Rm)
            meta["tau"] = np.mean(tau)
            meta["Cm"] = np.mean(Cm)
        ipsc_df.loc[i,'type'] = 'ipsc'
        ipsc_df.loc[i,"Ri"] = meta["Ri"]
        ipsc_df.loc[i,"Rm"] = meta["Rm"]
        ipsc_df.loc[i,"Rs"] = meta["Rs"]
        ipsc_df.loc[i,"tau"] = meta["tau"]  
        ipsc_df.loc[i,"Cm"] = meta["Cm"]
        
    os.chdir(epsc_loc) 
    for i,x in epsc_df.iterrows():
        Ri, Rs, Rm, tau, Cm = [], [], [],  [], []
        data, sr = read_abf(x.file_loc)
        for start in np.arange(0, len(data), sr):
            baseline = np.mean(data[start + baseline_start:start + baseline_end])
            peak     = np.min(data[start + peak_start:start + peak_end])
            inp      = np.mean(data[start + inp_start:start + inp_end])
            
            input_resistance = -0.01/((inp - baseline)*1e-12)
            series_resistance = -0.01/((peak - baseline)*1e-12)
            membrane_resistance = input_resistance - series_resistance 
            input_resistance = input_resistance / 1e6
            membrane_resistance = membrane_resistance / 1e6
            series_resistance = series_resistance / 1e6
            
            tau_perc = 1 - (1-1/np.exp(1))
            y_data = data[start+cm_peak_start:start + cm_base_start]
            x_data = np.arange(0,len(y_data),1)
            temp = np.max(y_data) - np.mean(data[start + int(0.9*sr):start+int(0.95*sr)])
            Itau = tau_perc * temp
            y_data2 = np.full(len(x_data),Itau)
            
            plt.plot(x_data,y_data,'-')
            plt.plot(x_data,y_data2,'-')
            idx = np.argwhere(np.diff(np.sign(y_data - y_data2))).flatten()
            plt.plot(x_data[idx],y_data[idx],'ro')
            
            try:
                Tau = idx[len(idx)-1] * 1000/sr       
            except:
                pass
            
            cm = (Tau / 1000)/ (membrane_resistance * 1e6) * 1e12
            
            Ri.append(input_resistance)
            Rs.append(series_resistance)
            Rm.append(membrane_resistance)
            tau.append(Tau)
            Cm.append(cm)
            meta["Ri"] = np.mean(Ri)
            meta["Rs"] = np.mean(Rs)
            meta["Rm"] = np.mean(Rm)
            meta["tau"] = np.mean(tau)
            meta["Cm"] = np.mean(Cm)
        epsc_df.loc[i,'type'] = 'epsc'
        epsc_df.loc[i,"Ri"] = meta["Ri"]
        epsc_df.loc[i,"Rm"] = meta["Rm"]
        epsc_df.loc[i,"Rs"] = meta["Rs"]
        epsc_df.loc[i,"tau"] = meta["tau"]  
        epsc_df.loc[i,"Cm"] = meta["Cm"]                 
    return ipsc_df, epsc_df

def export_data():
    ipsc_df, epsc_df = test()
    all_data = ipsc_df.append(epsc_df)
    return all_data

def seal_stats(excel_df):

    excel_df = excel_df.loc[excel_df.Cm > 20]
    vars        = ['Ri','Rs','Rm','Cm']
    for v in vars:
        print(v, tstd(excel_df.loc[excel_df.group == "pos", v]), 'Z+ standard deviation')
        print(v, tstd(excel_df.loc[excel_df.group == "neg", v]), 'Z- standard deviation')
        print(v, levene(excel_df.loc[excel_df.group == "pos", v], excel_df.loc[excel_df.group == "neg", v]))
        print(v, ttest_ind(excel_df.loc[excel_df.group == "pos", v], excel_df.loc[excel_df.group == "neg", v]), "\n")