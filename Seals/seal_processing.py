import numpy as np
import pandas as pd
import os
from neo.rawio import AxonRawIO
from numpy import trapz

### This script uses two resources to automatically analyse seal tests
### ABF files and a master excel spreadsheet is needed



#Function that reads ABF files and returns all datapoints and sampling rate
def read_abf(file_loc):
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header["signal_channels"][0][2])
    data = a.get_analogsignal_chunk() * a.header["signal_channels"][0][5]
    return data[:, 0], sr

#Reads a masterfile that contains information about cells, adapts cellname to full filename 
def read_masterfile():
    meta_loc = '//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Metadata'
    os.chdir(meta_loc)
    ipsc_df = pd.read_excel('mipsc_metafile.xlsx') 
    epsc_df = pd.read_excel('mepsc_metafile.xlsx')
    for i,x in ipsc_df.iterrows():
        file_loc = x['name'] + ' ' + str(x.cell) + '-' + str(x.seal)
        file_loc = file_loc + '.ABF'
        ipsc_df.loc[i,'file_loc'] = file_loc
        
        ipsc_df.loc[ipsc_df['mouse'].str.contains('F'), 'gender'] = 'Female'
        ipsc_df.loc[~ipsc_df['mouse'].str.contains('F'), 'gender'] = 'Male'
        
        epsc_df.loc[epsc_df['mouse'].str.contains('F'), 'gender'] = 'Female'
        epsc_df.loc[~epsc_df['mouse'].str.contains('F'), 'gender'] = 'Male'
            
    for i,x in epsc_df.iterrows():
        file_loc = x['name'] + ' ' + str(x.cell) + '-' + str(x.seal)
        file_loc = file_loc + '.ABF'
        epsc_df.loc[i,'file_loc'] = file_loc
    
    vars = ['II','III','IV-V'] #AZ
    for i,v in enumerate(vars):
        ipsc_df.loc[ipsc_df['lobule'] == v,'domain'] = 0
        epsc_df.loc[epsc_df['lobule'] == v,'domain'] = 0
    vars = ['VI','VII'] #CZ
    for i,v in enumerate(vars):
        ipsc_df.loc[ipsc_df['lobule'] == v,'domain'] = 1
        epsc_df.loc[epsc_df['lobule'] == v,'domain'] = 1
    vars = ['VIII','IX'] #PZ
    for i,v in enumerate(vars):
        ipsc_df.loc[ipsc_df['lobule'] == v,'domain'] = 2
        epsc_df.loc[epsc_df['lobule'] == v,'domain'] = 2
    
    
    #['Paramedian','Crus1','Crus2','Simple']
        
    return ipsc_df, epsc_df

#Uses the previous two functions
def obtain_data():
    ipsc_df, epsc_df = read_masterfile()
    ipsc_file_loc = master_sheet.file_loc 
    epsc_file_loc = master_sheet.file_loc
    
    ABF_loc = '//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/mIPSCs/Combined/seals/abf' 
    os.chdir(ABF_loc) 
    ipsc_data = pd.DataFrame()
    for i in ipsc_file_loc:
        data, sr = read_abf(i)
        print(data)
        data = pd.Series(data)
        ipsc_data = ipsc_data.append(data, ignore_index=True)
    epsc_data = pd.DataFrame()
    for o in epsc_file_loc:
        data, sr = read_abf(i)
        data = pd.Series(data)
        epsc_data = epsc_data.append(data, ignore_index=True)
    return ipsc_data, epsc_data


def analyse_data():
    ipsc_df, epsc_df = read_masterfile()
    sr = 50000
    
    ipsc_loc = '//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/mIPSCs/Combined/seals/abf'
    epsc_loc = '//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/mEPSCs/Combined/seals/abf'
 
    baseline_start   = int(0.65 * sr)
    baseline_end     = int(0.7 * sr)
    peak_start       = int(0.7 * sr)
    peak_end         = int(0.725 * sr)
    inp_start        = int(0.780 * sr)
    inp_end          = int(0.8 * sr)
    meta = {}
    
    os.chdir(ipsc_loc) 
    for i,x in ipsc_df.iterrows():
        data, sr = read_abf(x.file_loc)
        Ri, Rm, Rs, holding, Cm, Area = [],[],[],[],[],[] #Place where you store all values for every iteration
        for start in np.arange(0, len(data), sr):
            baseline = np.mean(data[start + baseline_start:start + baseline_end])
            peak     = np.min(data[start + peak_start:start + peak_end])
            inp      = np.mean(data[start + inp_start:start + inp_end])
            
            #resistance and holding
            input_resistance = -0.01/((inp - baseline)*1e-12)
            series_resistance = -0.01/((peak - baseline)*1e-12)
            membrane_resistance = input_resistance - series_resistance
            holding_current = baseline            
            input_resistance = input_resistance / 1e6
            membrane_resistance = membrane_resistance / 1e6
            series_resistance = series_resistance / 1e6
            
            #capacitance
            y_data = data[int(0.8*sr)+start:int(0.9 *sr)+start]
            y_data = y_data - np.min(y_data)
            y_data = y_data - np.min(y_data[4000:5000])
            x_data = np.arange(0,100,0.02)
            area = trapz(y_data, x_data)
            dv = 0.01 #V
            coulombs = area * 1e-15#C
            capacitance = (coulombs / dv) * 1e12 #
            Tau = (membrane_resistance * capacitance ) * 1e-3
            
            #This next part is important, since it will mean together multiple sweeps (so it does not matter how long your seal test is)
            Ri.append(input_resistance)
            Rs.append(series_resistance)
            Rm.append(membrane_resistance)
            holding.append(holding_current)
            Cm.append(capacitance)
            Area.append(area)
            meta["Ri"] = np.mean(Ri)
            meta["Rs"] = np.mean(Rs)
            meta["Rm"] = np.mean(Rm)
            meta["holding"] = holding_current * -1 #Remove/add the if you want inverted values
            meta["Cm"] = np.mean(Cm)
            meta["Area"] = np.mean(Area)
            meta["Tau"] = np.mean(Tau)
        ipsc_df.loc[i,'type'] = 'ipsc'
        ipsc_df.loc[i,"Ri"] = meta["Ri"]
        ipsc_df.loc[i,"Rm"] = meta["Rm"]
        ipsc_df.loc[i,"Rs"] = meta["Rs"]
        ipsc_df.loc[i,"holding"] = meta["holding"]
        ipsc_df.loc[i,'Cm'] = np.mean(Cm)
        ipsc_df.loc[i,'area'] = np.mean(Area)
        ipsc_df.loc[i,'Tau'] = meta['Tau']
    
    os.chdir(epsc_loc)
    for i,x in epsc_df.iterrows():
        data, sr = read_abf(x.file_loc)
        Ri, Rm, Rs, holding = [],[],[],[] #Place where you store all values for every iteration
        for start in np.arange(0, len(data), sr):
            baseline = np.mean(data[start + baseline_start:start + baseline_end])
            peak     = np.min(data[start + peak_start:start + peak_end])
            inp      = np.mean(data[start + inp_start:start + inp_end])
            
            #resistance and holding
            input_resistance = -0.01/((inp - baseline)*1e-12)
            series_resistance = -0.01/((peak - baseline)*1e-12)
            membrane_resistance = input_resistance - series_resistance
            holding_current = baseline            
            input_resistance = input_resistance / 1e6
            membrane_resistance = membrane_resistance / 1e6
            series_resistance = series_resistance / 1e6
                      
            #capacitance
            y_data = data[int(0.8*sr)+start:int(0.9 *sr)+start]
            y_data = y_data - np.min(y_data)
            y_data = y_data - np.min(y_data[4900:5000])
            x_data = np.arange(0,100,0.02)
            area = trapz(y_data, x_data)
            dv = 0.01 #V
            coulombs = area * 1e-15#C
            capacitance = (coulombs / dv) * 1e12 #pF
            Tau = (membrane_resistance * capacitance ) *1e-3
            
            #This next part is important, since it will mean together multiple sweeps (so it does not matter how long your seal test is)
            Ri.append(input_resistance)
            Rs.append(series_resistance)
            Rm.append(membrane_resistance)
            holding.append(holding_current)
            Cm.append(capacitance)
            Area.append(area)
            meta["Ri"] = np.mean(Ri)
            meta["Rs"] = np.mean(Rs)
            meta["Rm"] = np.mean(Rm)
            meta["holding"] = holding_current * -1 #Remove/add the if you want inverted values
            meta['Tau'] = np.mean(Tau)
        epsc_df.loc[i,'type'] = 'epsc'    
        epsc_df.loc[i,"Ri"] = meta["Ri"]
        epsc_df.loc[i,"Rm"] = meta["Rm"]
        epsc_df.loc[i,"Rs"] = meta["Rs"]
        epsc_df.loc[i,"holding"] = meta["holding"]
        epsc_df.loc[i,'Cm'] = np.mean(Cm)
        epsc_df.loc[i,'area'] = np.mean(Area)
        epsc_df.loc[i,'Tau'] = meta['Tau']
    all_data = ipsc_df.append(epsc_df)
    return ipsc_df, epsc_df, all_data

def exclusion():
    exclusion_df = pd.DataFrame()
    a,b,c = analyse_data()
    c.loc[:,'Ratio (Rs/Ri)'] = (c.Rs / c.Ri) *100
    
    #Criteria
    RsRi_cutoff = 15    
    Ri_cutoff    = 500
    holding_cutoff = 500
    
    #exclusion
    RsRi_df = c[c['Ratio (Rs/Ri)'] > RsRi_cutoff]
    Ri_df   = c[c['Ri'] > Ri_cutoff]
    holding_df = c[c['holding'] > holding_cutoff]
    
    exclusion_df = exclusion_df.append(RsRi_df)
    exclusion_df = exclusion_df.append(Ri_df)
    exclusion_df = exclusion_df.append(holding_df)
    exclusion_df = exclusion_df.drop_duplicates()
    
    index_to_drop = exclusion_df.index
    clean_data = c.drop(index_to_drop)
    
    return clean_data
    
def export_data():
    clean_data = exclusion()
    writer = pd.ExcelWriter('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Analysis_tools/output/seal_data_combined.xlsx')
    clean_data.to_excel(writer)
    writer.save()
    
    