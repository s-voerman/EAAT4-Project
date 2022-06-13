import numpy as np
import pandas as pd
import os
from neo.rawio import AxonRawIO
from numpy import trapz
import matplotlib.pyplot as plt
from scipy.stats import levene, ttest_ind
from scipy.stats import tstd
from scipy.stats import sem

def read_abf(file_loc):
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header["signal_channels"][0][2])
    data = a.get_analogsignal_chunk() * a.header["signal_channels"][0][5]
    return data[:, 0], sr

def read_masterfile():
    meta_loc = '//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Metadata'
    os.chdir(meta_loc)
    meta = pd.read_excel('test.xlsx') 

    for i,x in meta.iterrows():
        file_loc = x['name'] + ' ' + str(x.cell) + '-' + str(x.pre)
        file_loc = file_loc + '.ABF'
        meta.loc[i,'file_loc'] = file_loc

    return meta

def test():
    meta = read_masterfile()
    
    #ltd_loc = '//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/LTD/Data'
    ltd_loc = '//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/Test'
    
    sr=50000
    baseline_start   = int(0.65 * sr)
    baseline_end     = int(0.7 * sr)
    peak_start       = int(0.7 * sr)
    peak_end         = int(0.725 * sr)
    inp_start        = int(0.780 * sr)
    inp_end          = int(0.8 * sr)
    cm_peak_start    = int(0.8*sr)
    cm_base_start    = int(0.9 * sr)
    
    meta2 = {}
    output_df = pd.DataFrame()
    os.chdir(ltd_loc)
    for i,x in meta.iterrows():
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
            y_data = data[int(0.7*sr)+start:int(0.8 *sr)+start]
            y_data = y_data * - 1
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
            
            meta2["Ri"] = np.mean(Ri)
            meta2["Rs"] = np.mean(Rs)
            meta2["Rm"] = np.mean(Rm)
            meta2["holding"] = holding_current * -1 #Remove/add the if you want inverted values
            meta2["Cm"] = np.mean(Cm)
            meta2["Area"] = np.mean(Area)
            meta2['Tau'] = np.mean(Tau)
        meta.loc[i,"Ri"] = meta2["Ri"]
        meta.loc[i,"Rm"] = meta2["Rm"]
        meta.loc[i,"Rs"] = meta2["Rs"]
        meta.loc[i,"holding"] = meta2["holding"]
        meta.loc[i,'Cm'] = np.mean(Cm)
        meta.loc[i,'Tau'] = np.mean(Tau)
        meta.loc[i,'area'] = np.mean(Area)
    
    return meta

def export_data():
    ipsc_df, epsc_df = test()
    all_data = ipsc_df.append(epsc_df)
    return all_data

def export_data(data):
    output_df = pd.DataFrame()
    vars = ['Ri','Rs','Rm','Cm','Tau']
    for v in vars:
        pos_data = data.loc[data.group == 'pos',v]
        neg_data = data.loc[data.group == 'neg',v]
        
        SEM_neg = sem(neg_data)
        std_neg = np.std(neg_data)
        mean_neg = np.mean(neg_data)
        median_neg = np.median(neg_data)
        
        SEM_pos = sem(pos_data)
        std_pos = np.std(pos_data)
        mean_pos = np.mean(pos_data)
        median_pos = np.median(pos_data)
        
        output_df.loc['neg_mean',v] = mean_neg
        output_df.loc['neg_SEM',v] = SEM_neg
        output_df.loc['neg_median',v] = median_neg
        output_df.loc['neg_sd',v] = std_neg

        output_df.loc['pos_mean',v] = mean_pos
        output_df.loc['pos_SEM',v] = SEM_pos
        output_df.loc['pos_median',v] = median_pos
        output_df.loc['pos_sd',v] = std_pos
    
    output_df.to_excel(os.path.join(output_loc, "passivepropertiesPT.xlsx"))
    return output_df
    

def seal_stats(a):
    vars        = ['Ri','Rs','Rm','Cm','Tau']
    a = test()
    for v in vars:
        print(v, tstd(a.loc[a.group == "pos", v]), 'Z+ standard deviation')
        print(v, tstd(a.loc[a.group == "neg", v]), 'Z- standard deviation')
        print(v, levene(a.loc[a.group == "pos", v], a.loc[a.group == "neg", v]))
        print(v, ttest_ind(a.loc[a.group == "pos", v], a.loc[a.group == "neg", v]), "\n")
        
if __name__ == "__main__":
    base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
    output_loc = os.path.join(base_folder, "Analysis_tools/output/")