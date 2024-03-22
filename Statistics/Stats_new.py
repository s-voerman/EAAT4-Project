import os
import pandas as pd
import numpy as np
from scipy.stats import levene, ttest_ind, tstd, sem, shapiro, mannwhitneyu
from itertools import combinations

def plasticity_stats(data, output_loc):
    data = data.loc[data.time > -6]
    stat_df = pd.DataFrame()
    all_data = data.groupby(['group','cell','pre_post']).median().reset_index() 
    for i, x in all_data.iterrows():
        all_data['INDEX'] = all_data['pre_post'] + '_' + all_data['group']
    
    vars = ['norm_EPSC1_amplitude', 'norm_PPR', 'PPR','r_cap', 'norm_r_cap','r_membrane', 
            'norm_r_membrane', 'holding', 'norm_holding']
    groups = ['pre_pos','pre_neg','post_pos','post_neg']
    cc = list(combinations(groups,2))
    for i in range(len(cc)):
        a = all_data.loc[all_data.INDEX == cc[i][0]] #First group in pairing
        b = all_data.loc[all_data.INDEX == cc[i][1]] #Second group in pairing
        for v in vars:
            temp_df = pd.DataFrame()
            avg1, avg2, sem1, sem2, std1, std2, sw1, sw2 = [],[],[],[],[],[],[],[]
            tt, l , mw = [], [], []
            
            mean1 = np.mean(a.loc[:,v])
            mean2 = np.mean(b.loc[:,v])
            sderr1 = sem(a.loc[:,v])
            sderr2 = sem(b.loc[:,v])
            sd1 = tstd(a.loc[:,v])
            sd2 = tstd(b.loc[:,v])
            sh1 = shapiro(a.loc[:,v]) 
            sh2 = shapiro(b.loc[:,v])
            sh1 = sh1[1]
            sh2 = sh2[1]
            
            avg1.append(mean1)
            avg2.append(mean2)
            sem1.append(sderr1)
            sem2.append(sderr2)
            std1.append(sd1)
            std2.append(sd2)
            sw1.append(sh1)
            sw2.append(sh2)
            
            T = ttest_ind(a.loc[:,v],b.loc[:,v])
            MW = mannwhitneyu(a.loc[:,v],b.loc[:,v])
            L = levene(a.loc[:,v],b.loc[:,v])
            tt.append(T[1])
            l.append(L[1])
            mw.append(MW[1])

            temp_df.loc[:,'Mean1'] = pd.Series(avg1)
            temp_df.loc[:,'Mean2'] = pd.Series(avg2)
            temp_df.loc[:,'SEM1'] = pd.Series(sem1)
            temp_df.loc[:,'SEM2'] = pd.Series(sem2)
            temp_df.loc[:,'SD1'] = pd.Series(std1)
            temp_df.loc[:,'SD2'] = pd.Series(std2)
            temp_df.loc[:,'Sh-W1'] = pd.Series(sw1)
            temp_df.loc[:,'Sh-W2'] = pd.Series(sw2)               
            temp_df.loc[:,'TTest'] = pd.Series(tt)
            temp_df.loc[:,'Mann-Whitney'] = pd.Series(mw)
            temp_df.loc[:,'Levene'] = pd.Series(l)
            
            temp_df.loc[:,'Group1'] = str(cc[i][0]) 
            temp_df.loc[:,'Group2'] = str(cc[i][1])
            temp_df.loc[:,'Variable'] = str(v)
            temp_df.loc[:,'Temp'] = str(cc[i][0] + ' ' + str(cc[i][1]))
            stat_df = stat_df.append(temp_df)
            stat_df = stat_df.loc[stat_df.Temp != 'pre_neg post_pos']
            stat_df = stat_df.loc[stat_df.Temp != 'pre_pos post_neg']
            stat_df = stat_df.drop(columns =['Temp'])
    cols_to_order = ['Group1','Group2', 'Variable']
    new_columns = cols_to_order + (stat_df.columns.drop(cols_to_order).tolist())
    stat_df = stat_df[new_columns]
    
    stat_df.to_excel(os.path.join(output_loc, 'LTP_PT_statistics.xlsx'))
    return stat_df
    
def descriptives_only(data):
    data = data.loc[data.time >-6]
    all_data = data.groupby(['group','cell','pre_post']).median().reset_index() 
    
    stat_df = pd.DataFrame()
    vars = ['norm_EPSC1_amplitude', 'norm_PPR', 'PPR','r_cap','r_membrane', 'holding']
    groups = ['pos', 'neg']
    pre_post = ['pre', 'post']
    for z in groups:
        for p in pre_post: 
            y = all_data.loc[(all_data.group == z) & (all_data.pre_post == p)]
            for v in vars:
                temp_df = pd.DataFrame()
                avg, ser, std = [], [], []
                mean = np.mean(y.loc[:,v])
                sderr = sem(y.loc[:,v])
                sd = tstd(y.loc[:,v])
        
                avg.append(mean)
                ser.append(sderr)
                std.append(sd)
                temp_df.loc[:,'Mean'] = pd.Series(avg)
                temp_df.loc[:,'SEM'] = pd.Series(ser)
                temp_df.loc[:,'SD'] = pd.Series(std)
                temp_df['group'] = str(z)
                temp_df['pre_post'] = str(p)
                temp_df['variable'] = str(v)
                stat_df = stat_df.append(temp_df)

    return stat_df

def plasticity_stats_zones(data, output_loc):
    data = data.loc[(data.time > 19) & (data.time < 25)] #for this analysis I only need data from 20-25 mins
    data = data.groupby(['group','cell','domain']).median().reset_index()
    data.loc[(data.domain == 0), 'domain'] = 'AZ'
    data.loc[(data.domain == 1), 'domain'] = 'CZ'
    data.loc[(data.domain == 2), 'domain'] = 'PZ'
    stat_df = pd.DataFrame()
    
    vars = ['norm_EPSC1_amplitude']
    groups = ['pos', 'neg']
    domain_groups = ['pos_AZ', 'neg_AZ', 'pos_CZ', 'neg_CZ', 'pos_PZ', 'neg_PZ']
    cc = list(combinations(domain_groups,2))
    for i, x in data.iterrows():
        data['INDEX'] = data['group'] + '_' + data['domain']
    for i in range(len(cc)):
        a = data.loc[data.INDEX == cc[i][0]] #First group in pairing
        b = data.loc[data.INDEX == cc[i][1]] #Second group in pairing
        for v in vars:
            temp_df = pd.DataFrame()
            avg1, avg2, sem1, sem2, std1, std2, sw1, sw2 = [],[],[],[],[],[],[],[]
            tt, l , mw = [], [], []
        
            mean1 = np.mean(a.loc[:,v])
            mean2 = np.mean(b.loc[:,v])
            sderr1 = sem(a.loc[:,v])
            sderr2 = sem(b.loc[:,v])
            sd1 = tstd(a.loc[:,v])
            sd2 = tstd(b.loc[:,v])
            try:
                sh1 = shapiro(a.loc[:,v]) 
                sh2 = shapiro(b.loc[:,v])
                sh1 = sh1[1]
                sh2 = sh2[1]
            except:
                pass
            
            
            avg1.append(mean1)
            avg2.append(mean2)
            sem1.append(sderr1)
            sem2.append(sderr2)
            std1.append(sd1)
            std2.append(sd2)
            try:
                sw1.append(sh1)
                sw2.append(sh2)
            except:
                pass
            
            T = ttest_ind(a.loc[:,v],b.loc[:,v])
            MW = mannwhitneyu(a.loc[:,v],b.loc[:,v])
            L = levene(a.loc[:,v],b.loc[:,v])
            tt.append(T[1])
            l.append(L[1])
            mw.append(MW[1])

            temp_df.loc[:,'Mean1'] = pd.Series(avg1)
            temp_df.loc[:,'Mean2'] = pd.Series(avg2)
            temp_df.loc[:,'SEM1'] = pd.Series(sem1)
            temp_df.loc[:,'SEM2'] = pd.Series(sem2)
            temp_df.loc[:,'SD1'] = pd.Series(std1)
            temp_df.loc[:,'SD2'] = pd.Series(std2)
            temp_df.loc[:,'Sh-W1'] = pd.Series(sw1)
            temp_df.loc[:,'Sh-W2'] = pd.Series(sw2)
            temp_df.loc[:,'Mann-Whitney'] = pd.Series(tt)            
            temp_df.loc[:,'TTest'] = pd.Series(tt)
            temp_df.loc[:,'Levene'] = pd.Series(l)
            
            temp_df.loc[:,'Group1'] = str(cc[i][0]) 
            temp_df.loc[:,'Group2'] = str(cc[i][1])
            temp_df.loc[:,'Variable'] = str(v)
            stat_df = stat_df.append(temp_df)
    cols_to_order = ['Group1','Group2', 'Variable']
    new_columns = cols_to_order + (stat_df.columns.drop(cols_to_order).tolist())
    stat_df = stat_df[new_columns]
    
    for v in vars:
        pos_data = data.loc[data.group == groups[0]]
        neg_data = data.loc[data.group == groups[1]]
        temp_df = pd.DataFrame()
        avg1, avg2, sem1, sem2, std1, std2, sw1, sw2 = [],[],[],[],[],[],[],[]
        tt, l , mw = [], [], []
        
        mean1 = np.mean(pos_data.loc[:,v])
        mean2 = np.mean(neg_data.loc[:,v])
        sderr1 = sem(pos_data.loc[:,v])
        sderr2 = sem(neg_data.loc[:,v])
        sd1 = tstd(pos_data.loc[:,v])
        sd2 = tstd(neg_data.loc[:,v])
        try:
            sh1 = shapiro(pos_data.loc[:,v]) 
            sh2 = shapiro(neg_data.loc[:,v])
            sh1 = sh1[1]
            sh2 = sh2[1]
        except:
            pass
            
        
        avg1.append(mean1)
        avg2.append(mean2)
        sem1.append(sderr1)
        sem2.append(sderr2)
        std1.append(sd1)
        std2.append(sd2)
        try:
            sw1.append(sh1)
            sw2.append(sh2)
        except:
            pass
            
        T = ttest_ind(pos_data.loc[:,v],neg_data.loc[:,v])
        MW = mannwhitneyu(pos_data.loc[:,v],neg_data.loc[:,v])
        L = levene(pos_data.loc[:,v],neg_data.loc[:,v])
        tt.append(T[1])
        l.append(L[1])
        mw.append(MW[1])

        temp_df.loc[:,'Mean1'] = pd.Series(avg1)
        temp_df.loc[:,'Mean2'] = pd.Series(avg2)
        temp_df.loc[:,'SEM1'] = pd.Series(sem1)
        temp_df.loc[:,'SEM2'] = pd.Series(sem2)
        temp_df.loc[:,'SD1'] = pd.Series(std1)
        temp_df.loc[:,'SD2'] = pd.Series(std2)
        temp_df.loc[:,'Sh-W1'] = pd.Series(sw1)
        temp_df.loc[:,'Sh-W2'] = pd.Series(sw2)
        temp_df.loc[:,'Mann-Whitney'] = pd.Series(tt)            
        temp_df.loc[:,'TTest'] = pd.Series(tt)
        temp_df.loc[:,'Levene'] = pd.Series(l)
        
        temp_df.loc[:,'Group1'] = 'pos' 
        temp_df.loc[:,'Group2'] = 'neg'
        temp_df.loc[:,'Variable'] = str(v)
        stat_df = stat_df.append(temp_df)
    cols_to_order = ['Group1','Group2', 'Variable']
    new_columns = cols_to_order + (stat_df.columns.drop(cols_to_order).tolist())
    stat_df = stat_df[new_columns]
    
    stat_df.to_excel(os.path.join(output_loc, 'LTP_zones.xlsx'))        
    return stat_df

def epsc_stats(e, output_loc):
    stat_df = pd.DataFrame()
    groups = ['pos', 'neg']
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    for v in vars:
        a = e.loc[e.Zebrin == groups[0]] #Positive group
        b = e.loc[e.Zebrin == groups[1]] #Negative group
        temp_df = pd.DataFrame()
        avg1, avg2, sem1, sem2, std1, std2, sw1, sw2 = [],[],[],[],[],[],[],[]
        tt, l , mw = [], [], []
        
        mean1 = np.mean(a.loc[:,v])
        mean2 = np.mean(b.loc[:,v])
        sderr1 = sem(a.loc[:,v])
        sderr2 = sem(b.loc[:,v])
        sd1 = tstd(a.loc[:,v])
        sd2 = tstd(b.loc[:,v])
        sh1 = shapiro(a.loc[:,v]) 
        sh2 = shapiro(b.loc[:,v])
        sh1 = sh1[1]
        sh2 = sh2[1]
            
        avg1.append(mean1)
        avg2.append(mean2)
        sem1.append(sderr1)
        sem2.append(sderr2)
        std1.append(sd1)
        std2.append(sd2)
        sw1.append(sh1)
        sw2.append(sh2)
            
        T = ttest_ind(a.loc[:,v],b.loc[:,v])
        MW =  mannwhitneyu(a.loc[:,v],b.loc[:,v])
        L = levene(a.loc[:,v],b.loc[:,v])
        tt.append(T[1])
        l.append(L[1])
        mw.append(MW[1])
        
        temp_df.loc[:,'Mean1'] = pd.Series(avg1)
        temp_df.loc[:,'Mean2'] = pd.Series(avg2)
        temp_df.loc[:,'SEM1'] = pd.Series(sem1)
        temp_df.loc[:,'SEM2'] = pd.Series(sem2)
        temp_df.loc[:,'SD1'] = pd.Series(std1)
        temp_df.loc[:,'SD2'] = pd.Series(std2)   
        temp_df.loc[:,'Sh-W1'] = pd.Series(sw1)
        temp_df.loc[:,'Sh-W2'] = pd.Series(sw2)          
        temp_df.loc[:,'TTest'] = pd.Series(tt)
        temp_df.loc[:,'Mann-Whitney'] = pd.Series(mw)
        temp_df.loc[:,'Levene'] = pd.Series(l)
            
        temp_df.loc[:,'Group1'] = str(groups[0]) 
        temp_df.loc[:,'Group2'] = str(groups[1])
        temp_df.loc[:,'Variable'] = str(v)
        stat_df = stat_df.append(temp_df)
    cols_to_order = ['Group1','Group2', 'Variable']
    new_columns = cols_to_order + (stat_df.columns.drop(cols_to_order).tolist())
    stat_df = stat_df[new_columns]
    
    stat_df.to_excel(os.path.join(output_loc, 'epsc_stats.xlsx'))        
    return stat_df

def ipsc_stats(i, output_loc):
    stat_df = pd.DataFrame()
    groups = ['pos', 'neg']
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    for v in vars:
        a = i.loc[i.Zebrin == groups[0]] #Positive group
        b = i.loc[i.Zebrin == groups[1]] #Negative group
        temp_df = pd.DataFrame()
        avg1, avg2, sem1, sem2, std1, std2, sw1, sw2 = [],[],[],[],[],[],[],[]
        tt, l , mw = [], [], []
        
        mean1 = np.mean(a.loc[:,v])
        mean2 = np.mean(b.loc[:,v])
        sderr1 = sem(a.loc[:,v])
        sderr2 = sem(b.loc[:,v])
        sd1 = tstd(a.loc[:,v])
        sd2 = tstd(b.loc[:,v])
        sh1 = shapiro(a.loc[:,v]) 
        sh2 = shapiro(b.loc[:,v])
        sh1 = sh1[1]
        sh2 = sh2[1]
            
        avg1.append(mean1)
        avg2.append(mean2)
        sem1.append(sderr1)
        sem2.append(sderr2)
        std1.append(sd1)
        std2.append(sd2)
        sw1.append(sh1)
        sw2.append(sh2)
            
        T = ttest_ind(a.loc[:,v],b.loc[:,v])
        MW =  mannwhitneyu(a.loc[:,v],b.loc[:,v])
        L = levene(a.loc[:,v],b.loc[:,v])
        tt.append(T[1])
        l.append(L[1])
        mw.append(MW[1])
        
        temp_df.loc[:,'Mean1'] = pd.Series(avg1)
        temp_df.loc[:,'Mean2'] = pd.Series(avg2)
        temp_df.loc[:,'SEM1'] = pd.Series(sem1)
        temp_df.loc[:,'SEM2'] = pd.Series(sem2)
        temp_df.loc[:,'SD1'] = pd.Series(std1)
        temp_df.loc[:,'SD2'] = pd.Series(std2)   
        temp_df.loc[:,'Sh-W1'] = pd.Series(sw1)
        temp_df.loc[:,'Sh-W2'] = pd.Series(sw2)          
        temp_df.loc[:,'TTest'] = pd.Series(tt)
        temp_df.loc[:,'Mann-Whitney'] = pd.Series(mw)
        temp_df.loc[:,'Levene'] = pd.Series(l)
            
        temp_df.loc[:,'Group1'] = str(groups[0]) 
        temp_df.loc[:,'Group2'] = str(groups[1])
        temp_df.loc[:,'Variable'] = str(v)
        stat_df = stat_df.append(temp_df)
    cols_to_order = ['Group1','Group2', 'Variable']
    new_columns = cols_to_order + (stat_df.columns.drop(cols_to_order).tolist())
    stat_df = stat_df[new_columns]
    
    stat_df.to_excel(os.path.join(output_loc, 'ipsc_stats.xlsx'))        
    return stat_df

def epsc_stats_zones(e, output_loc):
    data = e.groupby(['Zebrin','Cellname','domain']).median().reset_index()
    data.loc[(data.domain == 0), 'domain'] = 'AZ'
    data.loc[(data.domain == 1), 'domain'] = 'CZ'
    data.loc[(data.domain == 2), 'domain'] = 'PZ'
    for j, x in data.iterrows():
        data['INDEX'] = data['Zebrin'] + '_' + data['domain']
        
    stat_df = pd.DataFrame()
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    domain_groups = ['pos_AZ', 'neg_AZ', 'pos_CZ', 'neg_CZ', 'pos_PZ', 'neg_PZ']
    cc = list(combinations(domain_groups,2))
    

    for j in range(len(cc)):
        a = data.loc[data.INDEX == cc[j][0]] #First group in pairing
        b = data.loc[data.INDEX == cc[j][1]] #Second group in pairing
        for v in vars:
            temp_df = pd.DataFrame()
            avg1, avg2, sem1, sem2, std1, std2, sw1, sw2 = [],[],[],[],[],[],[],[]
            tt, l , mw = [], [], []
        
            mean1 = np.mean(a.loc[:,v])
            mean2 = np.mean(b.loc[:,v])
            sderr1 = sem(a.loc[:,v])
            sderr2 = sem(b.loc[:,v])
            sd1 = tstd(a.loc[:,v])
            sd2 = tstd(b.loc[:,v])
            try:
                sh1 = shapiro(a.loc[:,v]) 
                sh2 = shapiro(b.loc[:,v])
            except:
                pass
            try:
                sh1 = sh1[1]
                sh2 = sh2[1]
            except:
                pass
        
            avg1.append(mean1)
            avg2.append(mean2)
            sem1.append(sderr1)
            sem2.append(sderr2)
            std1.append(sd1)
            std2.append(sd2)
            sw1.append(sh1)
            sw2.append(sh2)
            
            T = ttest_ind(a.loc[:,v],b.loc[:,v])
            MW =  mannwhitneyu(a.loc[:,v],b.loc[:,v])
            L = levene(a.loc[:,v],b.loc[:,v])
            tt.append(T[1])
            l.append(L[1])
            mw.append(MW[1])
            
            temp_df.loc[:,'Mean1'] = pd.Series(avg1)
            temp_df.loc[:,'Mean2'] = pd.Series(avg2)
            temp_df.loc[:,'SEM1'] = pd.Series(sem1)
            temp_df.loc[:,'SEM2'] = pd.Series(sem2)
            temp_df.loc[:,'SD1'] = pd.Series(std1)
            temp_df.loc[:,'SD2'] = pd.Series(std2)   
            temp_df.loc[:,'Sh-W1'] = pd.Series(sw1)
            temp_df.loc[:,'Sh-W2'] = pd.Series(sw2)          
            temp_df.loc[:,'TTest'] = pd.Series(tt)
            temp_df.loc[:,'Mann-Whitney'] = pd.Series(mw)
            temp_df.loc[:,'Levene'] = pd.Series(l)
            
            temp_df.loc[:,'Group1'] = str(cc[j][0]) 
            temp_df.loc[:,'Group2'] = str(cc[j][1])
            temp_df.loc[:,'Variable'] = str(v)
            stat_df = stat_df.append(temp_df)
    cols_to_order = ['Group1','Group2', 'Variable']
    new_columns = cols_to_order + (stat_df.columns.drop(cols_to_order).tolist())
    stat_df = stat_df[new_columns]
    stat_df.to_excel(os.path.join(output_loc, 'epsc_zones.xlsx'))   
    return stat_df

def ipsc_stats_zones(i, output_loc):
    data = i.groupby(['Zebrin','Cellname','domain']).median().reset_index()
    data.loc[(data.domain == 0), 'domain'] = 'AZ'
    data.loc[(data.domain == 1), 'domain'] = 'CZ'
    data.loc[(data.domain == 2), 'domain'] = 'PZ'
    for j, x in data.iterrows():
        data['INDEX'] = data['Zebrin'] + '_' + data['domain']
        
    stat_df = pd.DataFrame()
    vars = ["amplitude", "Frequency (Hz)", "rise", "tau"]
    domain_groups = ['pos_AZ', 'neg_AZ', 'pos_CZ', 'neg_CZ', 'pos_PZ', 'neg_PZ']
    cc = list(combinations(domain_groups,2))
    

    for j in range(len(cc)):
        a = data.loc[data.INDEX == cc[j][0]] #First group in pairing
        b = data.loc[data.INDEX == cc[j][1]] #Second group in pairing
        for v in vars:
            temp_df = pd.DataFrame()
            avg1, avg2, sem1, sem2, std1, std2, sw1, sw2 = [],[],[],[],[],[],[],[]
            tt, l , mw = [], [], []
        
            mean1 = np.mean(a.loc[:,v])
            mean2 = np.mean(b.loc[:,v])
            sderr1 = sem(a.loc[:,v])
            sderr2 = sem(b.loc[:,v])
            sd1 = tstd(a.loc[:,v])
            sd2 = tstd(b.loc[:,v])
            try:
                sh1 = shapiro(a.loc[:,v]) 
                sh2 = shapiro(b.loc[:,v])
            except:
                pass
            try:
                sh1 = sh1[1]
                sh2 = sh2[1]
            except:
                pass
        
            avg1.append(mean1)
            avg2.append(mean2)
            sem1.append(sderr1)
            sem2.append(sderr2)
            std1.append(sd1)
            std2.append(sd2)
            sw1.append(sh1)
            sw2.append(sh2)
            
            T = ttest_ind(a.loc[:,v],b.loc[:,v])
            MW =  mannwhitneyu(a.loc[:,v],b.loc[:,v])
            L = levene(a.loc[:,v],b.loc[:,v])
            tt.append(T[1])
            l.append(L[1])
            mw.append(MW[1])
            
            temp_df.loc[:,'Mean1'] = pd.Series(avg1)
            temp_df.loc[:,'Mean2'] = pd.Series(avg2)
            temp_df.loc[:,'SEM1'] = pd.Series(sem1)
            temp_df.loc[:,'SEM2'] = pd.Series(sem2)
            temp_df.loc[:,'SD1'] = pd.Series(std1)
            temp_df.loc[:,'SD2'] = pd.Series(std2)   
            temp_df.loc[:,'Sh-W1'] = pd.Series(sw1)
            temp_df.loc[:,'Sh-W2'] = pd.Series(sw2)          
            temp_df.loc[:,'TTest'] = pd.Series(tt)
            temp_df.loc[:,'Mann-Whitney'] = pd.Series(mw)
            temp_df.loc[:,'Levene'] = pd.Series(l)
            
            temp_df.loc[:,'Group1'] = str(cc[j][0]) 
            temp_df.loc[:,'Group2'] = str(cc[j][1])
            temp_df.loc[:,'Variable'] = str(v)
            stat_df = stat_df.append(temp_df)
    cols_to_order = ['Group1','Group2', 'Variable']
    new_columns = cols_to_order + (stat_df.columns.drop(cols_to_order).tolist())
    stat_df = stat_df[new_columns]
    stat_df.to_excel(os.path.join(output_loc, 'ipsc_zones.xlsx'))   
    return stat_df

if __name__ == "__main__":
    base_folder = ""
    data_loc = os.path.join(base_folder, "Data/LTP/Combined/abf")
    meta_loc = os.path.join(base_folder, "Metadata/ltpmeta_select_strict - Copy.xlsx")
    output_loc = os.path.join(base_folder, "Analysis_tools/output/")
