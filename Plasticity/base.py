import matplotlib as mpl
import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numpy import trapz
import helper_functions as hf
import pandas as pd
from scipy.stats import levene, ttest_ind, tstd, sem, pearsonr

plt.ion()
pd.options.display.max_columns = 20
pd.options.display.width = 160
pd.options.display.max_rows = 130

#Function for reading supplemental data.
#Make sure the excel sheets are named according to the code below, or change the code to reflect the names of the excel sheets.
#Make sure the excel sheets are in the correct folder.
def read_supplements():
    ltd = pd.read_excel('ltd_data.xlsx')
    ltd_pairs = pd.read_excel('ltd_pairs.xlsx')
    ltp = pd.read_excel('ltp_data.xlsx')
    ltp_pt = pd.read_excel('ltp_data_pt.xlsx')
    return ltd, ltd_pairs, ltp, ltp_pt

class process_single_cell(object):
    def __init__(self, file_loc, output_loc, cell, group, pre, ind, post, row, lobule, confirmed, domain):
        self.file_loc = file_loc
        self.meta_row = row
        self.output_loc = output_loc
        self.cell = cell
        self.group = group
        self.pre = pre
        self.ind = ind
        self.post = post
        self.lobule = lobule
        self.confirmed = confirmed
        self.domain = domain
        # timepoints
        self.epsc1 = 0.1
        self.epsc2 = 0.2
        self.seal = [0.7, 0.8]
        self.sweep_size = 1.0
        return

    def process_file(self, data, sr, prepost):
        output = {}
        sweep_size = int(sr * self.sweep_size)
        windows = np.arange(0, len(data), sweep_size)
        output_pd = pd.DataFrame()
        print(len(data), sr, len(data)/sr)
        for i, x in enumerate(windows):
            output = {}
            d = data[x:x+sweep_size]
            output["pre_post"] = prepost
            if prepost == "post":
                output["time"] = i + 15
            else:
                output["time"] = i
            for epsc, start in zip(["EPSC1", "EPSC2"], [int(self.epsc1*sr),
                                                        int(self.epsc2*sr)]):
                baseline = np.mean(d[start-400:start-10])
                peak = np.min(d[start+1:start+400])
                amplitude = peak - baseline
                for sv, v in zip(["baseline", "peak", "amplitude"],
                                 [baseline, peak, amplitude]):
                    output["_".join([epsc, sv])] = v
            seal_baseline = np.mean(d[int(0.6*sr):int(0.695*sr)])
            seal_cap = np.min(d[int(0.69*sr):int(0.71*sr)]) - seal_baseline
            seal_input = np.mean(d[int(0.77*sr):int(0.795*sr)]) - seal_baseline
            output["r_cap"] = (-10e-3 / (seal_cap*1e-12)) / 1e6
            output["r_input"] = (-10e-3 / (seal_input*1e-12)) / 1e6
            output["r_membrane"] = output["r_input"] - output["r_cap"]
            output["holding"] = baseline
            output_pd = output_pd.append(output, ignore_index=True)
        output_pd.loc[:, "cell"] = self.cell
        output_pd.loc[:, "group"] = self.group
        output_pd.loc[:, "lobule"] = self.lobule
        output_pd.loc[:, "confirmed"] = self.confirmed
        output_pd.loc[:, "domain"] = self.domain
        return output_pd
    
    def complex_spikes(self, data, sr):
        output = []
        window = np.arange(0, len(data), sr)
        for start in window:
            # baseline = np.min(data[int(0.3 * sr)+start:int(0.4 * sr)+start])
            # cs_data = data[int(0.08 * sr)+start:int(0.12 * sr)+start]
            baseline = np.min(data[int(0.3 * sr)+start:int(0.4 * sr)+start])
            cs_data = data[int(0.20 * sr)+start:int(0.24 * sr)+start]
            cs_data = cs_data- baseline
            x_data = np.arange(0,40,0.02)
            area = trapz(cs_data, x_data)
            output.append(area)
        cs_data = np.mean(output)
        return cs_data
    
    def run(self):
        pre_data = np.array([])
        for f in self.pre:
            filename = self.cell + "-" + f + ".ABF"
            data, sr = hf.read_abf(os.path.join(self.file_loc, filename))
            pre_data = np.concatenate([pre_data, data])
        post_data = np.array([])
        for f in self.post:
            filename = self.cell + "-" + f + ".ABF"
            data, sr = hf.read_abf(os.path.join(self.file_loc, filename))
            post_data = np.concatenate([post_data, data])
        cs_data = np.array([])
        for f in self.ind:
            filename = self.cell + "-" + f + ".ABF"
            data, sr = hf.read_abf(os.path.join(self.file_loc, filename))
            cs_data = np.concatenate([cs_data, data])    
        sr = 50000
        pre_data = self.process_file(pre_data, sr, "pre")
        pre_data.loc[:, "time"] = [x - (np.max(pre_data.time)) for x in pre_data.time]
        pre_data.loc[:, "time"] = [int(x/3)-1 for x in pre_data.time]
        post_data = self.process_file(post_data, sr, "post")
        post_data.loc[:, "time"] = [int(x/3) for x in post_data.time]
        cs_data = self.complex_spikes(cs_data, sr)
        output_data = pd.concat([pre_data, post_data])
        output_data['CS'] = cs_data
        #normalize
        #I still need to make a choice about whether to normalize from -10 to 0 or from -5 to 0. 
        output_data.loc[:, "PPR"] = output_data.EPSC2_amplitude / output_data.EPSC1_amplitude
        for i, v in enumerate(["EPSC1_amplitude", "EPSC2_amplitude",
                               "r_cap", "r_membrane", "PPR", "holding"]):
            output_data.loc[:, "norm_"+v] = (output_data.loc[:, v] /
                                            np.mean(output_data.loc[(output_data.time > -6) & (output_data.time < 1), v]))
            #print(np.mean(output_data.loc[(output_data.time > -10) & (output_data.time < 0), v]))
        return output_data
    

def find_meta(meta_loc, data_loc):
    all_data = pd.DataFrame()
    meta = pd.read_excel(meta_loc)
    
    vars = ['II','III','IV'] #AZ
    for i,v in enumerate(vars):
        meta.loc[meta['lobule'] == v,'domain'] = 0
    vars = ['V','VI'] #CZ
    for i,v in enumerate(vars):
        meta.loc[meta['lobule'] == v,'domain'] = 1
    vars = ['VII','VIII','IX','X'] #PZ
    for i,v in enumerate(vars):
        meta.loc[meta['lobule'] == v,'domain'] = 2
        
    for i, row in meta.iterrows():
        lobule = row.lobule
        confirmed = row.confirmed
        domain = row.domain
        group = row.group
        cell = " ".join([str(row["name"]), str(row["cell"])])
        file_prefix = row.cell
        if isinstance(row.pre, str):
            pre = row.pre.split(", ")
        else:
            pre = [str(row.pre)]
        if isinstance(row.post, str):
            post = row.post.split(", ")
        else:
            post = [str(row.post)]
        if isinstance(row.ind, str):
            ind = row.ind.split(", ")
        else:
            ind = [str(row.ind)]
        p = process_single_cell(data_loc, "./output", cell, group, pre, ind, post, row, lobule, confirmed, domain)
        output = p.run()

        all_data = pd.concat([all_data, output])
    return all_data


def overview_graph(data, output_loc):
    #data = pd.DataFrame(data.groupby(["cell", "group", "time"]).mean().reset_index())
    plt.style.use('seaborn-whitegrid')
    data = pd.DataFrame(data.groupby(["cell", "time"]).mean().reset_index())
    for cell in data.cell.unique():
        d = data.loc[data.cell == cell]
        fig, ax = plt.subplots(6,2, sharex=True, figsize=(6,12))
        fig.suptitle(cell)
        labels = {"EPSC1_amplitude" : "EPSC1",
                  "EPSC2_amplitude" : "EPSC2",
                  "r_cap" : "Rseries",
                  "r_membrane" : "Rmembrane",
                  "PPR": "PPR",
                  "holding":"holding"}
        for i, v in enumerate(["EPSC1_amplitude", "EPSC2_amplitude",
                               "r_cap", "r_membrane", "PPR", "holding"]):
            ax[i,0].scatter(d.time, d[v], s=4)
            ax[i,0].set_xlabel("Time")
            ax[i,0].set_ylabel(labels[v])

            ax[i,1].scatter(d.time, d["norm_"+v], s=4)
            ax[i,1].axhline(0.75, lw=0.7,ls='--')
            ax[i,1].axhline(1.25, lw=0.7,ls='--')
            ax[i,1].axvline(25, lw=0.3, color="black")
            ax[i,1].axvline(-5, lw=0.3, color="black")
            ax[i,1].axvline(0, lw=0.3, color="black")
            ax[i,1].axvline(30, lw=0.3, color="black")
            ax[i,1].set_ylim(0,2)
            ax[i,1].set_xlabel("Time")
            ax[i,1].set_ylabel("norm" + labels[v])

        plt.tight_layout()
        plt.subplots_adjust(top=.92)
        fig.savefig(os.path.join(output_loc, "summary_"+cell+".png"))
        plt.close()
    return

def ltp_stats(data, output_loc):
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
    stat_df.to_excel(os.path.join(output_loc, 'LTD_stats.xlsx'))
    return stat_df


def LTP_plot(data):
    #sns.set(style = 'ticks')
    data = data.loc[data.norm_PPR < 5]
    cm = 1/2.54
    pre_neg = data.loc[(data.group == "neg") & (data.time < 0) & (data.time > -6)]
    pre_pos = data.loc[(data.group == "pos") & (data.time < 0) & (data.time > -6)]
    pre_neg.loc[:,'time'] = pre_neg.time + 0.5
    
    post_neg = data.loc[(data.group == "neg") & (data.time >0) & (data.time < 36)]
    post_pos = data.loc[(data.group == "pos") & (data.time >0) & (data.time < 36)]
    post_neg.loc[:,'time'] = post_neg.time + 0.5
    fig, ax = plt.subplots(1,1, figsize=(0.75*(18.3*cm),3))
    pre_data = pre_neg.append(pre_pos)
    post_data = post_neg.append(post_pos)
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=pre_data,ax=ax, err_style="bars", ci=95,style='group',
                 markers=['o','^'],hue_order=['pos','neg'], legend=False, linewidth=0.5, markersize=3,err_kws={'elinewidth':0.5},
                 dashes=False)
  
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=post_data,ax=ax, err_style="bars", ci=95, style='group',
                 markers=['o','^'],hue_order=['pos','neg'], legend=False, linewidth=0.5, markersize=3,err_kws={'elinewidth':0.5},
                 dashes=False)
    
    ax.axhline(y=1.0,xmin=-10.0,xmax=35.5, linestyle = '--', color='black', lw=0.5)
    ax.axvline(x=0, linestyle='--', lw=0.5, color='black')
    ax.axvline(x=5, linestyle='--', lw=0.5, color='black')
    ax.set_ylabel('EPSC (%)')
    ax.set_xlabel('Time (minutes)')
    ax.set_xlim([-6,37])
    ax.set_ylim([0.6,1.4])
    ax.set_xticks([-5,0,5,10,15,20,25,30,35])
    ax.set_yticks([0.6,0.8,1.0,1.2,1.4])
    ax.set_yticklabels([60,80,100,120,140])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines["bottom"].set_color('black')
    ax.spines["left"].set_color('black')
    ax.set_facecolor("white")
    #ax.legend([],[],frameon=False)
    ax.tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
    ax.tick_params(axis='y',which='both',left=True,right=False,labelbottom=True)
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'LTD_final3.png'), dpi=1000)
    return fig, ax


def flocculus_plot(data):
    cm = 1/2.54
    pre_floc = data.loc[(data.group == "flocculus") & (data.time < 0) & (data.time > -11)]
    post_floc = data.loc[(data.group == "flocculus") & (data.time > 0) & (data.time < 36)]
    pre_floc.loc[:,'time'] = pre_floc.time + 0.25
    post_floc.loc[:,'time'] = post_floc.time + 0.25
    pre_neg = data.loc[(data.group == "neg") & (data.time < 0) & (data.time > -11)]
    pre_pos = data.loc[(data.group == "pos") & (data.time < 0) & (data.time > -11)]
    pre_neg.loc[:,'time'] = pre_neg.time + 0.5
    post_neg = data.loc[(data.group == "neg") & (data.time >0) & (data.time < 36)]
    post_pos = data.loc[(data.group == "pos") & (data.time >0) & (data.time < 36)]
    post_neg.loc[:,'time'] = post_neg.time + 0.5
    pre_data = pre_neg.append(pre_pos)
    pre_data = pre_data.append(pre_floc)
    post_data = post_neg.append(post_pos)
    post_data = post_data.append(post_floc)
    
    fig, ax = plt.subplots(1,1, figsize=(18.3*cm,4), dpi = 300)
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=pre_data,ax=ax, err_style="bars", ci=95,style='group',
                 markers=['o','^','X'],hue_order=['pos','neg','flocculus'], 
                 err_kws={'capsize':5,'elinewidth':2,'capthick':2})
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=post_data,ax=ax, err_style="bars", ci=95, style='group',
                 markers=['o','^','X'],hue_order=['pos','neg','flocculus'], legend=False, 
                 err_kws={'capsize':5,'elinewidth':2,'capthick':2})
    ax.axhline(y=1.0,xmin=-10.0,xmax=35.5, linestyle = '--', color='black', lw=0.5)
    ax.set_ylabel('Normalized EPSC1 Amplitude')
    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel('EPSC (%)')
    ax.set_xlabel('Time (minutes)')
    ax.set_xlim([-12,37])
    ax.set_ylim([0,1.8])
    ax.set_xticks([-10,-5,0,5,10,15,20,25,30,35])
    ax.set_yticklabels([0,20,40,60,80,100,120,140,160,180])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines["bottom"].set_color('black')
    ax.spines["left"].set_color('black')
    ax.set_facecolor("white")
    ax.tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
    ax.tick_params(axis='y',which='both',left=True,right=False,labelbottom=True)
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'ltd_plot_flocculus.png'),dpi = 1000)
    return fig, ax

def test_plot(data):
    data = data.loc[data.norm_PPR < 5] #removes one very unlikely outlier that ruins the figure 
    cm = 1/2.54
    pre_neg = data.loc[(data.group == "neg") & (data.time < 0) & (data.time > -11)]
    pre_pos = data.loc[(data.group == "pos") & (data.time < 0) & (data.time > -11)]
    pre_neg.loc[:,'time'] = pre_neg.time + 0.5
    
    post_neg = data.loc[(data.group == "neg") & (data.time >0) & (data.time < 36)]
    post_pos = data.loc[(data.group == "pos") & (data.time >0) & (data.time < 36)]
    post_neg.loc[:,'time'] = post_neg.time + 0.5
    fig, ax = plt.subplots(1,2, figsize=(0.7*(18.3*cm),2.5))
    ax[0].axhline(y=1.0,xmin=-10.0,xmax=35.5, linestyle = '--', color='black', lw=0.5)
    ax[1].axhline(y=1.0,xmin=-10.0,xmax=35.5, linestyle = '--', color='black', lw=0.5)
    ax[0].axvline(x=0, linestyle='--', lw=0.5, color='black')
    ax[0].axvline(x=5, linestyle='--', lw=0.5, color='black')
    ax[1].axvline(x=0, linestyle='--', lw=0.5, color='black')
    ax[1].axvline(x=5, linestyle='--', lw=0.5, color='black')
    pre_data = pre_neg.append(pre_pos)
    post_data = post_neg.append(post_pos)
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=pre_data,ax=ax[0], err_style="bars", ci=95,style='group', 
                 markers=['o','^'], markersize=3, linewidth=0.5, hue_order=['pos','neg'], dashes=False,
                 err_kws={'elinewidth':0.5})
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=post_data,ax=ax[0], err_style="bars", ci=95, style='group',
                 markers=['o','^'], markersize=3, linewidth=0.5,hue_order=['pos','neg'], legend=False, 
                 err_kws={'elinewidth':0.5}, dashes=False)
    sns.lineplot(x='time',y='norm_PPR', hue='group', data=pre_data,ax=ax[1], err_style="bars", ci=95,style='group',
                 hue_order=['pos','neg'], markers=['o','^'], markersize=3, linewidth=0.5, err_kws={'elinewidth':0.5},dashes=False)
    sns.lineplot(x='time',y='norm_PPR', hue='group', data=post_data,ax=ax[1], err_style="bars", ci=95, style='group',
                 markers=['o','^'], markersize=3, dashes=False,
                 linewidth=0.5, legend=False, hue_order=['pos','neg'], err_kws={'elinewidth':0.5})    
    ax[0].set_ylabel('EPSC (%)')
    ax[0].set_xlabel('Time (minutes)')
    ax[1].set_ylabel('PPR (%)')
    ax[1].set_xlabel('Time (minutes)')
    ax[0].set_xlim([-6,37])
    ax[0].set_ylim([0.6,1.4])
    ax[0].set_xticks([-5,0,5,10,15,20,25,30,35])
    ax[0].set_yticks([0.6,0.8,1.0,1.2,1.4,1.6])
    ax[0].set_yticklabels([60,80,100,120,140,160])
    ax[1].set_xlim([-6,37])
    ax[1].set_ylim([0.6,1.4])
    ax[1].set_xticks([-5,0,5,10,15,20,25,30,35])
    ax[1].set_yticks([0.6,0.8,1.0,1.2,1.4,1.6])
    ax[1].set_yticklabels([60,80,100,120,140,160])
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines["bottom"].set_color('black')
    ax[0].spines["left"].set_color('black')
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines["bottom"].set_color('black')
    ax[1].spines["left"].set_color('black')
    ax[0].legend().set_visible(False)
    ax[1].legend().set_visible(False)
    ax[0].set_facecolor("white")
    ax[1].set_facecolor("white")
    ax[0].tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
    ax[0].tick_params(axis='y',which='both',left=True,right=False,labelbottom=True)
    ax[1].tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
    ax[1].tick_params(axis='y',which='both',left=True,right=False,labelbottom=True)
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'LTD_final.pdf'))
    return fig, ax

def cp_plot(data):
    cm=1/2.54
    pre_data = data.loc[data.pre_post == 'pre']
    post_data = data.loc[data.pre_post == 'post']
    pre_data = pre_data.groupby(['cell','group']).mean().reset_index()
    post_data = post_data.groupby(['cell','group']).mean().reset_index()
    all_data = data.groupby(['pre_post','cell','group']).mean().reset_index()
    all_data['holding'] = all_data['holding'] * -1
    
    fig, ax = plt.subplots(2,4, figsize = ((13.6*cm),13.6*cm))
    vars = ['PPR','norm_PPR','r_cap','norm_r_cap','r_membrane','norm_r_membrane','holding','norm_holding']
    for i, v in enumerate(vars):
        ax_x = int(i/4)
        if i < 4:
            ax_y = i
        else:
            ax_y = i-4
        sns.boxplot(x='pre_post',y=v,hue='group',data=all_data,ax=ax[ax_x,ax_y],order=['pre','post'], hue_order=['pos','neg'],
                fliersize=0, linewidth=0.5)
        sns.stripplot(x='pre_post',y=v,hue='group',data=all_data,ax=ax[ax_x,ax_y],dodge=True,color='white', edgecolor='black', 
                      linewidth=0.5,
                      hue_order=['pos','neg'], order =['pre','post'], size=3, marker='^', jitter=0)
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.legend().remove()
        a.set_xlabel(" ")
        a.set_xticklabels(['Pre','Post'])
    ax[0,0].set_ylim(1.0,2.0)
    ax[0,1].set_ylim(0.4,1.6)
    ax[0,1].set_ylabel('Normalized PPR')
    ax[0,2].set_ylim(0,20)
    ax[0,2].set_ylabel('Series Resistance (MΩ)')
    ax[0,3].set_ylim(0.4,1.6)
    ax[0,3].set_ylabel('Normalized Series')
    ax[1,0].set_ylim(0,200)
    ax[1,0].set_ylabel('Membrane Resistance (MΩ)')    
    ax[1,1].set_ylim(0.4,1.6)
    ax[1,1].set_ylabel('Normalized Membrane')
    ax[1,2].set_ylim(0,600)
    ax[1,2].set_ylabel('Holding Current (pA)')
    ax[1,3].set_ylim(0.4,1.6)
    ax[1,3].set_ylabel('Normalized Holding')
    fig.tight_layout()
    fig.align_ylabels()
    fig.savefig(os.path.join(output_loc, 'misc_plasticity_info.pdf'))
    return fig, ax

def domain_plot(data):
    cm = 1/2.54
    neg_data = data.loc[(data.domain == 2) & (data.group == 'neg')]
    pos_data = data.loc[(data.domain == 0) & (data.group == 'pos')]
    plot_data = neg_data.append(pos_data)
    pre_data = plot_data.loc[(plot_data.time < 1) & (plot_data.time > -11)]
    post_data = plot_data.loc[(plot_data.time >  1) & (plot_data.time < 36)]
    fig, ax = plt.subplots(1,1, figsize=(18.3*cm,4), dpi = 300)
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=pre_data,ax=ax, err_style="bars", ci=95,style='group',
                 markers=['o','^'],hue_order=['pos','neg'], err_kws={'capsize':5,'elinewidth':2,'capthick':2})
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=post_data,ax=ax, err_style="bars", ci=95, style='group',
                 markers=['o','^'],hue_order=['pos','neg'], legend=False, err_kws={'capsize':5,'elinewidth':2,'capthick':2})
    ax.set_ylabel('Normalized EPSC1 Amplitude')
    ax.set_xlabel('Time (minutes)')
    ax.axhline(y=1.0,xmin=-10.0,xmax=35.5, linestyle = '--', color='black')
    ax.axvline(x=0.0,ymin=0,ymax=1, linestyle='--',color='black')
    ax.axvline(x=5.0,ymin=0,ymax=1, linestyle='--',color='black')
    ax.set_xlim([-12,37])
    ax.set_xticks([-10,-5,0,5,10,15,20,25,30,35])
    ax.set_yticks([0.4,0.6,0.8,1.0,1.2,1.4])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines["bottom"].set_color('black')
    ax.spines["left"].set_color('black')
    ax.set_facecolor("white")
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'regionmatch_plot.png'),dpi = 1000)
    return fig, ax

def extended_data_fig(data, pair_data, output_loc):
    all_data = data.loc[(data.time > 19) & (data.time < 26)] #Data for general graph, not necessary to have loc
    all_data = all_data.groupby(['group','cell','pre_post']).mean().reset_index()
    pair_data = pair_data.loc[(pair_data.time > 19) & (pair_data.time <26)]
    pair_data = pair_data.groupby(['group','cell','pre_post','lobule']).mean().reset_index()
    pair_data = pair_data.sort_values(by='group', ascending=False)
    palette1 = sns.color_palette('husl', 5, 0.7)
    cm=1/2.54
    fig, ax = plt.subplots(1,2, dpi=300, figsize=(0.3*(18.3*cm),3))
    sns.boxplot(data=all_data, x='group', y='norm_EPSC1_amplitude', fliersize=0, ax=ax[0], order=['pos','neg'], linewidth=0.5)
    sns.stripplot(data=all_data, x='group', y='norm_EPSC1_amplitude', ax=ax[0], order=['pos','neg'], size=4, hue='group',
                  hue_order=['pos','neg'], jitter=0.05, marker='^', edgecolor='black', color='white', linewidth=0.5)
    
    sns.lineplot(data=pair_data, x='group', y='norm_EPSC1_amplitude', hue='lobule', ax=ax[1], ci =None, palette=palette1,
                 sort=False)
    sns.swarmplot(data=pair_data, x='group', y='norm_EPSC1_amplitude', hue='lobule', ax=ax[1], edgecolor='black', linewidth=0.5,
                  palette=palette1, order=['pos','neg'])
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.set_ylim(0.4, 1.2)
        a.set_yticks([0.4,0.6,0.8,1.0,1.2])
        a.set_yticklabels([40,60,80,100,120])
    ax[0].set_ylabel('EPSC (%) (t=20-25 min)')
    ax[1].set_ylabel(' ')
    ax[0].legend().set_visible(False)
    ax[1].legend().set_visible(False)
    ax[0].set_xticklabels(['+','-'])
    ax[1].set_xticklabels(['+','-'])
    ax[0].set_xlabel('EAAT4')
    ax[1].set_xlabel('EAAT4')
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "extendedLTDdata.pdf"), dpi = 300)
    return fig, ax

def paired_plot(data, PT_data):
    data = data.loc[(data.time > 19) & (data.time < 26)] #Data for general graph, not necessary to have loc
    data = data.groupby(['group','cell','pre_post']).mean().reset_index()
    PT_data = PT_data.loc[(PT_data.time > 19) & (PT_data.time < 26)] #Data for general graph, not necessary to have loc
    PT_data = PT_data.groupby(['group','cell','pre_post']).mean().reset_index()
    cm = 1/2.54
    fig, ax = plt.subplots(1,2, figsize = (0.5*(18.3*cm), 3))
    sns.boxplot(data=data, x='group', y='norm_EPSC1_amplitude', fliersize=0, ax=ax[0], order=['pos','neg'], linewidth=0.5)
    sns.stripplot(data=data, x='group', y='norm_EPSC1_amplitude', ax=ax[0], order=['pos','neg'], size=4, hue='group',
                  hue_order=['pos','neg'], jitter=0.05, marker='^', edgecolor='black', color='white', linewidth=0.5)
    sns.boxplot(data=PT_data, x='group', y='norm_EPSC1_amplitude', fliersize=0, ax=ax[1], order=['pos','neg'], linewidth=0.5)
    sns.stripplot(data=PT_data, x='group', y='norm_EPSC1_amplitude', ax=ax[1], order=['pos','neg'], size=4, hue='group',
                  hue_order=['pos','neg'], jitter=0.05, marker='^', edgecolor='black', color='white', linewidth=0.5)
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.set_ylim(0.4, 2.0)
        a.set_yticks([0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0])
        a.set_yticklabels([40,60,80,100,120,140,160,180,200])
    ax[0].set_ylabel('EPSC (%) (t=20-25 min)')
    ax[1].set_ylabel(' ')
    ax[0].legend().set_visible(False)
    ax[1].legend().set_visible(False)
    ax[0].set_xticklabels(['+','-'])
    ax[1].set_xticklabels(['+','-'])
    ax[0].set_xlabel('EAAT4')
    ax[1].set_xlabel('EAAT4')
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "extended_LTP_data.pdf"))
    return

def seperate_plot(data):
    data = data.loc[data.norm_PPR < 5]
    cm = 1/2.54
    pre_neg = data.loc[(data.group == "neg") & (data.time < 0) & (data.time > -6)]
    pre_pos = data.loc[(data.group == "pos") & (data.time < 0) & (data.time > -6)]
    pre_neg.loc[:,'time'] = pre_neg.time + 0.5
    post_neg = data.loc[(data.group == "neg") & (data.time >0) & (data.time < 36)]
    post_pos = data.loc[(data.group == "pos") & (data.time >0) & (data.time < 36)]
    post_neg.loc[:,'time'] = post_neg.time + 0.5

    fig, ax = plt.subplots(1,2, figsize=(0.7*(18.3*cm),3))
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', data=pre_pos,ax=ax[0], err_style="bars", ci=95,style='group',
                 legend=False, linewidth=0.5, markersize=3,err_kws={'elinewidth':0.5}, markers=['^'],
                 dashes=False, color="#5BE12D")
    
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', data=post_pos,ax=ax[0], err_style="bars", ci=95,style='group',
                 legend=False, linewidth=0.5, markersize=3,err_kws={'elinewidth':0.5}, markers=['^'],
                 dashes=False, color="#5BE12D")

    sns.lineplot(x='time',y='norm_EPSC1_amplitude', data=pre_neg,ax=ax[1], err_style="bars", ci=95, style='group',
                 legend=False, linewidth=0.5, markersize=3,err_kws={'elinewidth':0.5}, markers=['o'],
                 dashes=False, color="#EC008C")
    
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', data=post_neg,ax=ax[1], err_style="bars", ci=95, style='group',
                 legend=False, linewidth=0.5, markersize=3,err_kws={'elinewidth':0.5}, markers=['o'],
                 dashes=False, color="#EC008C")
    
    for a in ax.flatten():
        a.axhline(y=1.0,xmin=-10.0,xmax=35.5, linestyle = '--', color='black', lw=0.5)
        a.axvline(x=0, linestyle='--', lw=0.5, color='black')
        a.axvline(x=5, linestyle='--', lw=0.5, color='black')
        a.set_ylabel('EPSC (%)')
        a.set_xlabel('Time (minutes)')
        a.set_xlim([-6,37])
        a.set_ylim([0.6,1.4])
        a.set_xticks([-5,0,5,10,15,20,25,30,35])
        a.set_yticks([0.6,0.8,1.0,1.2,1.4])
        a.set_yticklabels([60,80,100,120,140])
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.spines["bottom"].set_color('black')
        a.spines["left"].set_color('black')
        a.set_facecolor("white")
        #ax.legend([],[],frameon=False)
        a.tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
        a.tick_params(axis='y',which='both',left=True,right=False,labelbottom=True)
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'LTD_seperate.pdf'))
    return fig, ax

def undelayed_plot(data):
    all_data = data.loc[(data.time > 19) & (data.time < 26)] #Data for general graph, not necessary to have loc
    all_data = all_data.groupby(['group','cell','pre_post']).mean().reset_index()
    
    pre_neg = data.loc[(data.group == "neg") & (data.time < 0) & (data.time > -6)]
    pre_pos = data.loc[(data.group == "pos") & (data.time < 0) & (data.time > -6)]
    pre_neg.loc[:,'time'] = pre_neg.time + 0.5
    post_neg = data.loc[(data.group == "neg") & (data.time >0) & (data.time < 36)]
    post_pos = data.loc[(data.group == "pos") & (data.time >0) & (data.time < 36)]
    post_neg.loc[:,'time'] = post_neg.time + 0.5
    pre_data = pre_neg.append(pre_pos)
    post_data = post_neg.append(post_pos)
    cm=1/2.54
    fig, ax = plt.subplots(dpi=300, figsize=(0.3*(18.3*cm),3))
    sns.barplot(x='group',y='norm_EPSC1_amplitude', data=all_data, errorbar=('ci',95), width=0.5, errcolor='black', errwidth=0.5, capsize=0.3,
                facecolor=(0,0,0,0), edgecolor=(["#5BE12D","#EC008C"]), lw=0.5)
    sns.stripplot(data=all_data, x='group', y='norm_EPSC1_amplitude', ax=ax, order=['pos','neg'], size=4, hue='group',
                  hue_order=['pos','neg'], jitter=0, marker='^', edgecolor='black', color='white', linewidth=0.5)
    ax.axhline(y=1.0, linestyle='--', lw=0.5, color='black')
    ax.set_ylabel('EPSC (%)')
    ax.set_xlabel('EAAT4')
    ax.set_xticklabels(['+','-'])
    ax.set_yticks([0.4,0.6,0.8,1.0,1.2,1.4,1.6])
    ax.set_yticklabels([40,60,80,100,120,140,160])
    ax.set_ylim([0.4,1.4])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines["bottom"].set_color('black')
    ax.spines["left"].set_color('black')
    ax.legend().set_visible(False)
    ax.set_facecolor("white")
    ax.tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
    ax.tick_params(axis='y',which='both',left=True,right=False,labelbottom=True)
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'LTD_distribution.pdf'))
    return fig, ax

def export_data(data):
    #for RMAnova
    post_data = data.loc[data.pre_post == 'post']
    post_data = post_data.loc[post_data.time < 36]
    a = np.arange(4,30,5)
    for i in a:
        df = post_data.loc[(post_data.time > i) & (post_data.time <i+6)]
        #print(df)
        print(i+6)
        df = df.groupby(['cell','group']).mean().reset_index()
        df.to_excel('data' + str(i)+'.xlsx')
    return 

def cs_correlation(data):
    post_data = data.loc[data.pre_post == 'post']
    plot_data = post_data.groupby(['cell','group']).mean().reset_index()
    pos_data = plot_data.loc[plot_data.group == 'pos']
    neg_data = plot_data.loc[plot_data.group == 'neg']
    mean_data = post_data.groupby(['group']).mean().reset_index()
    fig, ax = plt.subplots()
    # sns.regplot(data = pos_data, x='CS', y='norm_EPSC1_amplitude',ax=ax, color='#5BE12D',
    #              scatter_kws={"s": 40, "color": 'black'},
    #              robust=False, ci=None, scatter=False)
    # sns.regplot(data = neg_data, x='CS', y='norm_EPSC1_amplitude',ax=ax, color='#EC008C',
    #             scatter_kws={"s": 40, "color": 'black'},
    #             robust=False, ci=None, scatter=False)
    #sns.lineplot(data = mean_data, x='CS', y='norm_EPSC1_amplitude', hue='group', color='black', ax=ax)
    #sns.scatterplot(data = plot_data, x='CS', y='norm_EPSC1_amplitude', hue='group',ax=ax, markers=['o','^'],style='group')
    ax.set_yticks([0.6,0.8,1.0,1.2])
    ax.set_yticklabels([60,80,100,120])
    ax.set_ylim([0.6,1.2])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('EPSC (%) - Post LTD induction')
    ax.set_xlabel('Complex spike charge')
    ax.legend().set_visible(False)
    fig.savefig(os.path.join(output_loc, 'LTD_cs_correlation.pdf'))
    return

def cs_correlation2(data):
    cm = 1/2.54
    post_data = data.loc[data.pre_post == 'post']
    plot_data = post_data.groupby(['cell','group']).mean().reset_index()
    pos_data = plot_data.loc[plot_data.group == 'pos']
    neg_data = plot_data.loc[plot_data.group == 'neg']
    pos_mean = np.mean(pos_data)
    pos_std = pos_data.sem()
    neg_mean = np.mean(neg_data)
    neg_std = neg_data.sem()
    fig, ax = plt.subplots(1,2, figsize=(0.75*(13.6*cm),3))
    sns.scatterplot(data = plot_data, x='CS', y='norm_EPSC1_amplitude', hue='group',ax=ax[0], markers=['^','o'],style='group', hue_order=['pos','neg'])
    ax[0].errorbar(x=pos_mean.CS, y=pos_mean.norm_EPSC1_amplitude, xerr=pos_std.CS, yerr=pos_std.norm_EPSC1_amplitude, elinewidth=2.5,
                capsize=5, capthick=2.5)
    ax[0].errorbar(x=neg_mean.CS, y=neg_mean.norm_EPSC1_amplitude, xerr=neg_std.CS, yerr=neg_std.norm_EPSC1_amplitude, elinewidth=2.5,
                capsize=5, capthick=2.5)
    
    data = data.loc[data.norm_PPR < 5]
    pre_neg = data.loc[(data.group == "neg") & (data.time < 0) & (data.time > -6)]
    pre_pos = data.loc[(data.group == "pos") & (data.time < 0) & (data.time > -6)]
    pre_neg.loc[:,'time'] = pre_neg.time + 0.5
    post_neg = data.loc[(data.group == "neg") & (data.time >0) & (data.time < 36)]
    post_pos = data.loc[(data.group == "pos") & (data.time >0) & (data.time < 36)]
    post_neg.loc[:,'time'] = post_neg.time + 0.5
    pre_data = pre_neg.append(pre_pos)
    post_data = post_neg.append(post_pos)
    sns.lineplot(x='time',y='norm_PPR', hue='group', data=pre_data,ax=ax[1], err_style="bars", ci=95,style='group',
                 markers=['o','^'],hue_order=['pos','neg'], legend=False, linewidth=0.5, markersize=3,err_kws={'elinewidth':0.5},
                 dashes=False)
  
    sns.lineplot(x='time',y='norm_PPR', hue='group', data=post_data,ax=ax[1], err_style="bars", ci=95, style='group',
                 markers=['o','^'],hue_order=['pos','neg'], legend=False, linewidth=0.5, markersize=3,err_kws={'elinewidth':0.5},
                 dashes=False)

    ax[0].set_yticks([0.6,0.8,1.0,1.2])
    ax[0].set_yticklabels([60,80,100,120])
    ax[0].set_ylim([0.6,1.2])
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].set_ylabel('EPSC (%) - Post LTD induction')
    ax[0].set_xlabel('Complex spike charge')
    ax[0].legend().set_visible(False)
    ax[1].axhline(y=1.0,xmin=-10.0,xmax=35.5, linestyle = '--', color='black', lw=0.5)
    ax[1].axvline(x=0, linestyle='--', lw=0.5, color='black')
    ax[1].axvline(x=5, linestyle='--', lw=0.5, color='black')
    ax[1].set_ylabel('PPR')
    ax[1].set_xlabel('Time (minutes)')
    ax[1].set_xlim([-6,37])
    ax[1].set_ylim([0.6,1.4])
    ax[1].set_xticks([-5,0,5,10,15,20,25,30,35])
    ax[1].set_yticks([0.6,0.8,1.0,1.2,1.4])
    ax[1].set_yticklabels([60,80,100,120,140])
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines["bottom"].set_color('black')
    ax[1].spines["left"].set_color('black')
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'LTD_PPR.pdf'))
    return

def cs_stats(data):
    post_data = data.loc[data.pre_post == 'post']
    plot_data = post_data.groupby(['cell','group']).mean().reset_index()
    pos_data = plot_data.loc[plot_data.group == 'pos']
    neg_data = plot_data.loc[plot_data.group == 'neg']
    print(levene(pos_data.norm_EPSC1_amplitude, neg_data.norm_EPSC1_amplitude), 'epsc')
    print(levene(pos_data.CS, neg_data.CS), 'CS')
    
    print(pearsonr(plot_data.norm_EPSC1_amplitude, plot_data.CS), 'all')
    print(pearsonr(pos_data.norm_EPSC1_amplitude, pos_data.CS), 'pos')
    print(pearsonr(neg_data.norm_EPSC1_amplitude, neg_data.CS), 'neg')
    return

def extract(data):
    data.to_excel(os.path.join(output_loc,'excel.xlsx'))
    return

if __name__ == "__main__":
    base_folder = ""
    data_loc = os.path.join(base_folder, "Data/LTD/Data")
    meta_loc = os.path.join(base_folder, "Metadata/ltdmeta_select_strict.xlsx")
    output_loc = os.path.join(base_folder, "Analysis_tools/output/")
    plt.style.use('seaborn-paper')
    #color = ["#FFFFFF","#808080"] #white/grey
    color = ["#5BE12D","#EC008C"] #Original green/red
    #color = ["#4DAF4A","#E41A1C"] #New green/red
    sns.set_palette(sns.color_palette(color), desat=1)
    #sns.set_palette('Set1',8,0.7)
    #sns.set_palette('Dark2',8)
    plt.rcParams['font.sans-serif'] = 'Arial'
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['xtick.minor.width'] = 0.3
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['ytick.minor.width'] = 0.3
    fonts = {"font.size":7, "axes.labelsize":7, "ytick.labelsize":7, "xtick.labelsize":7}
    plt.rcParams.update(fonts)
    # output = find_meta(meta_loc, data_loc)

