import pandas as pd
import numpy as np
import os

#This is mainly to prep excel or csv sheets for use in JASP

def convert_plasticity(data, output_loc):
    data = data.loc[data.time > -6]
    data = data.groupby(['group','cell','pre_post']).median().reset_index() 
    
    output = data[['cell','group','pre_post','norm_EPSC1_amplitude', 'norm_PPR','PPR','r_cap','r_membrane','holding']]
    output.to_excel(os.path.join(output_loc, 'LTP_PT_bayesian.xlsx'))
    output.to_csv(os.path.join(output_loc, 'LTP_PT_bayesian.csv'))
    return output

def convert_plasticity_time(data, output_loc):
    output = pd.DataFrame()
    post_data = data.loc[data.pre_post == 'post']
    post_data = post_data.loc[post_data.time < 36]
    
    unit = np.arange(4,30,5)
    for i in unit:
        df = post_data.loc[(post_data.time > i) & (post_data.time <i+6)]
        df = df.groupby(['cell','group']).median().reset_index()
        
        output['cell'] = pd.Series(df.cell)
        output['group'] = pd.Series(df.group)
        output[str(i + 1) + '-' + str(i + 6)] = pd.Series(df.norm_EPSC1_amplitude)
        
    output.to_excel(os.path.join(output_loc, 'LTP_PT_25_bayesian.xlsx'))
    output.to_csv(os.path.join(output_loc, 'LTP_PT_25_bayesian.csv'))        
    return output

def convert_plasticity_zones(data, output_loc):
    #NOTE: All data at time 20-25 mins
    output = pd.DataFrame()
    data = data.loc[(data.time > 19) & (data.time < 25)]
    data.loc[(data.domain == 0), 'domain'] = 'AZ'
    data.loc[(data.domain == 1), 'domain'] = 'CZ'
    data.loc[(data.domain == 2), 'domain'] = 'PZ'
    data = data.groupby(['group','cell','domain']).median().reset_index()
    
    
    groups = ['pos', 'neg']
    zones = ['AZ', 'CZ', 'PZ']
    for g in groups:
        for d in zones:
            temp = data.loc[(data.group == g) & (data.domain == d)]
            temp['domain_group'] = str(d) + '_' + str(g)
            output = output.append(temp)
            output = output[['cell','domain_group','norm_EPSC1_amplitude']]
            
    output.to_excel(os.path.join(output_loc, 'LTP_PT_zones_bayesian.xlsx'))
    output.to_csv(os.path.join(output_loc, 'LTP_PT_zones_bayesian.csv'))
    return output

def convert_psc(e, i, output_loc):
    vars = ['Cellname','Zebrin','amplitude', 'Frequency (Hz)','rise', 'tau']
    epsc_output = e[vars]
    ipsc_output = i[vars]
    
    epsc_output.to_excel(os.path.join(output_loc, 'epsc_bayesian.xlsx'))
    epsc_output.to_csv(os.path.join(output_loc, 'epsc_bayesian.csv'))
    ipsc_output.to_excel(os.path.join(output_loc, 'ipsc_bayesian.xlsx'))
    ipsc_output.to_csv(os.path.join(output_loc, 'ipsc_bayesian.csv'))
    return epsc_output, ipsc_output

def convert_psc_zones(e, i , output_loc):
    e.loc[(e.domain == 0), 'domain'] = 'AZ'
    e.loc[(e.domain == 1), 'domain'] = 'CZ'
    e.loc[(e.domain == 2), 'domain'] = 'PZ'
    i.loc[(i.domain == 0), 'domain'] = 'AZ'
    i.loc[(i.domain == 1), 'domain'] = 'CZ'
    i.loc[(i.domain == 2), 'domain'] = 'PZ'
    
    epsc_output = pd.DataFrame()
    ipsc_output = pd.DataFrame()
    
    vars = ['Cellname','domain_group','amplitude', 'Frequency (Hz)','rise', 'tau']
    groups = ['pos', 'neg']
    zones = ['AZ', 'CZ', 'PZ']
    for g in groups:
        for d in zones:
            temp = e.loc[(e.Zebrin == g) & (e.domain == d)]
            temp['domain_group'] = str(d) + '_' + str(g)
            epsc_output = epsc_output.append(temp)
            epsc_output = epsc_output[vars]
    for g in groups:
        for d in zones:
            temp = i.loc[(i.Zebrin == g) & (i.domain == d)]
            temp['domain_group'] = str(d) + '_' + str(g)
            ipsc_output = ipsc_output.append(temp)
            ipsc_output = ipsc_output[vars]
    
    epsc_output.to_excel(os.path.join(output_loc, 'epsc_zones_bayesian.xlsx'))
    epsc_output.to_csv(os.path.join(output_loc, 'epsc_zones_bayesian.csv'))
    ipsc_output.to_excel(os.path.join(output_loc, 'ipsc_zones_bayesian.xlsx'))
    ipsc_output.to_csv(os.path.join(output_loc, 'ipsc_zones_bayesian.csv'))
    return epsc_output, ipsc_output
    
    
    
        