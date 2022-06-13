import matplotlib as mpl
import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import helper_functions as hf
import pandas as pd
from scipy.stats import levene, ttest_ind
from scipy.stats import linregress
from scipy.stats import tstd
from scipy.stats import spearmanr
plt.ion()
pd.options.display.max_columns = 20
pd.options.display.width = 160
pd.options.display.max_rows = 130

class process_single_cell(object):
    def __init__(self, file_loc, output_loc, cell, group, pre, post, row, lobule, confirmed, domain):
        self.file_loc = file_loc
        self.meta_row = row
        self.output_loc = output_loc
        self.cell = cell
        self.group = group
        self.pre = pre
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

        sr = 50000
        pre_data = self.process_file(pre_data, sr, "pre")
        pre_data.loc[:, "time"] = [x - (np.max(pre_data.time)) for x in pre_data.time]
        pre_data.loc[:, "time"] = [int(x/3)-1 for x in pre_data.time]
        post_data = self.process_file(post_data, sr, "post")
        post_data.loc[:, "time"] = [int(x/3) for x in post_data.time]
        output_data = pd.concat([pre_data, post_data])
        
        #normalize
        output_data.loc[:, "PPR"] = output_data.EPSC2_amplitude / output_data.EPSC1_amplitude
        #output_data = output_data.loc[output_data.EPSC1_amplitude > -700] #temp
        for i, v in enumerate(["EPSC1_amplitude", "EPSC2_amplitude",
                               "r_cap", "r_membrane", "PPR", "holding"]):
            output_data.loc[:, "norm_"+v] = (output_data.loc[:, v] /
                                            np.mean(output_data.loc[(output_data.time > -10) & (output_data.time < 0), v]))
            #print(np.mean(output_data.loc[(output_data.time > -10) & (output_data.time < 0), v]))
        return output_data
    
def norm_test(data):
    pre_data = data.loc[data.pre_post == 'pre']
    post_data = data.loc[data.pre_post == 'post']
    
    pre_data = pre_data.groupby(['cell','time','pre_post','group']).mean().reset_index()
    post_data = post_data.groupby(['cell','time','pre_post','group']).mean().reset_index()
    output_data = pd.concat([pre_data, post_data])
    for i, v in enumerate(["EPSC1_amplitude", "EPSC2_amplitude",
                               "r_cap", "r_membrane", "PPR", "holding"]):
            output_data.loc[:, "norm2_"+v] = (output_data.loc[:, v] /
                                            np.mean(output_data.loc[(output_data.time > -10) & (output_data.time < 0), v]))
            
    pre_data2 = output_data.loc[output_data.pre_post == 'pre']
    post_data2 = output_data.loc[output_data.pre_post == 'post']
    fig, ax = plt.subplots()
    sns.lineplot(x='time',y='norm2_EPSC1_amplitude', hue='group', data=pre_data2,ax=ax, err_style="bars", ci=95,style='group',
                 markers=['o','^'],hue_order=['pos','neg'], err_kws={'capsize':5,'elinewidth':2,'capthick':2})
    sns.lineplot(x='time',y='norm2_EPSC1_amplitude', hue='group', data=post_data2,ax=ax, err_style="bars", ci=95, style='group',
                 markers=['o','^'],hue_order=['pos','neg'], legend=False, err_kws={'capsize':5,'elinewidth':2,'capthick':2})
    
    
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
        p = process_single_cell(data_loc, "./output", cell, group, pre, post, row, lobule, confirmed, domain)
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
            ax[i,1].axhline(0.8, lw=0.7,ls='--')
            ax[i,1].axhline(1.2, lw=0.7,ls='--')
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


def ltp_stats(data):
    plot_data  = data.groupby(['group','cell','pre_post']).median().reset_index()
    temp_data  = plot_data.groupby('pre_post')
    pre_data   = temp_data.get_group('pre')
    post_data  = temp_data.get_group('post')
    
    vars = ['norm_EPSC1_amplitude','norm_EPSC2_amplitude','PPR']
    for v in vars:
        print(v, tstd(pre_data.loc[pre_data.group == "pos", v]), 'Z+ standard deviation')
        print(v, tstd(pre_data.loc[pre_data.group == "neg", v]), 'Z- standard deviation')
        print(v, levene(pre_data.loc[pre_data.group == "pos", v], pre_data.loc[pre_data.group == "neg", v]))
        print(v, ttest_ind(pre_data.loc[pre_data.group == "pos", v], pre_data.loc[pre_data.group == "neg", v]), "\n")
        
    for v in vars:
        print(v, tstd(post_data.loc[post_data.group == "pos", v]), 'Z+ standard deviation')
        print(v, tstd(post_data.loc[post_data.group == "neg", v]), 'Z- standard deviation')
        print(v, levene(post_data.loc[post_data.group == "pos", v], post_data.loc[post_data.group == "neg", v]))
        print(v, ttest_ind(post_data.loc[post_data.group == "pos", v], post_data.loc[post_data.group == "neg", v]), "\n")
        
    return


def LTP_plot(data):
    #sns.set(style = 'ticks')
    cm = 1/2.54
    pre_neg = data.loc[(data.group == "neg") & (data.time < 0) & (data.time > -6)]
    pre_pos = data.loc[(data.group == "pos") & (data.time < 0) & (data.time > -6)]
    pre_neg.loc[:,'time'] = pre_neg.time + 0.5
    
    post_neg = data.loc[(data.group == "neg") & (data.time >0) & (data.time < 36)]
    post_pos = data.loc[(data.group == "pos") & (data.time >0) & (data.time < 36)]
    post_neg.loc[:,'time'] = post_neg.time + 0.5
    fig, ax = plt.subplots(1,1, figsize=(18.3*cm,4), dpi = 300)
    
    pre_data = pre_neg.append(pre_pos)
    post_data = post_neg.append(post_pos)


    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=pre_data,ax=ax, err_style="bars", ci=95,style='group',
                 markers=['o','^'],hue_order=['pos','neg'], err_kws={'capsize':5,'elinewidth':2,'capthick':2})
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=post_data,ax=ax, err_style="bars", ci=95, style='group',
                 markers=['o','^'],hue_order=['pos','neg'], legend=False, err_kws={'capsize':5,'elinewidth':2,'capthick':2})
    

    ax.axhline(y=1.0,xmin=-10.0,xmax=35.5, linestyle = '--', color='black', lw=0.5)

    
    ax.set_ylabel('Normalized EPSC1 Amplitude')
    ax.set_xlabel('Time (minutes)')
    
    ax.set_xlim([-12,37])
    ax.set_ylim([0.6,1.2])
    ax.set_xticks([-10,-5,0,5,10,15,20,25,30,35])
    ax.set_yticklabels([60,80,100,120])
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines["bottom"].set_color('black')
    ax.spines["left"].set_color('black')
    ax.set_facecolor("white")
    #ax.legend([],[],frameon=False)

    ax.tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
    ax.tick_params(axis='y',which='both',left=True,right=False,labelbottom=True)
    
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'ltd_plot.png'),dpi = 1000)
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
    fig, ax = plt.subplots(1,1, figsize=(18.3*cm,4), dpi = 300)
    
    pre_data = pre_neg.append(pre_pos)
    pre_data = pre_data.append(pre_floc)
    post_data = post_neg.append(post_pos)
    post_data = post_data.append(post_floc)
    
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
    data = data.loc[data.norm_PPR < 5] #removes one very unlikely datapoint 
    cm = 1/2.54
    pre_neg = data.loc[(data.group == "neg") & (data.time < 0) & (data.time > -11)]
    pre_pos = data.loc[(data.group == "pos") & (data.time < 0) & (data.time > -11)]
    pre_neg.loc[:,'time'] = pre_neg.time + 0.5
    
    post_neg = data.loc[(data.group == "neg") & (data.time >0) & (data.time < 36)]
    post_pos = data.loc[(data.group == "pos") & (data.time >0) & (data.time < 36)]
    post_neg.loc[:,'time'] = post_neg.time + 0.5
    fig, ax = plt.subplots(1,2, figsize=((0.7*18.3)*cm,3), dpi=300)
    
    ax[0].axhline(y=1.0,xmin=-10.0,xmax=35.5, linestyle = '--', color='black', lw=0.5)
    ax[1].axhline(y=1.0,xmin=-10.0,xmax=35.5, linestyle = '--', color='black', lw=0.5)
    ax[0].axvline(x=0, linestyle='--', lw=0.5, color='black')
    ax[0].axvline(x=5, linestyle='--', lw=0.5, color='black')
    ax[1].axvline(x=0, linestyle='--', lw=0.5, color='black')
    ax[1].axvline(x=5, linestyle='--', lw=0.5, color='black')
    
    pre_data = pre_neg.append(pre_pos)
    post_data = post_neg.append(post_pos)


    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=pre_data,ax=ax[0], err_style="bars", ci=95,style='group', 
                 markers=['o','^'], markersize=3, linewidth=0.5, hue_order=['pos','neg'],
                 err_kws={'capsize':2,'elinewidth':0.5,'capthick':0.5})
    sns.lineplot(x='time',y='norm_EPSC1_amplitude', hue='group', data=post_data,ax=ax[0], err_style="bars", ci=95, style='group',
                 markers=['o','^'], markersize=3, linewidth=0.5,hue_order=['pos','neg'], legend=False, 
                 err_kws={'capsize':2,'elinewidth':0.5,'capthick':0.5})
    sns.lineplot(x='time',y='norm_PPR', hue='group', data=pre_data,ax=ax[1], err_style="bars", ci=95,style='group',
                 hue_order=['pos','neg'], markers=['o','^'], markersize=3, linewidth=0.5, err_kws={'capsize':2,'elinewidth':0.5,'capthick':0.5})
    sns.lineplot(x='time',y='norm_PPR', hue='group', data=post_data,ax=ax[1], err_style="bars", ci=95, style='group',
                 markers=['o','^'], markersize=3,
                 linewidth=0.5, legend=False, hue_order=['pos','neg'], err_kws={'capsize':2,'elinewidth':0.5,'capthick':0.5})


    
    ax[0].set_ylabel('EPSC (%)')
    ax[0].set_xlabel('Time (minutes)')
    ax[1].set_ylabel('PPR (%)')
    ax[1].set_xlabel('Time (minutes)')
    
    ax[0].set_xlim([-12,37])
    ax[0].set_ylim([0.6,1.4])
    ax[0].set_xticks([-10,-5,0,5,10,15,20,25,30,35])
    ax[0].set_yticks([0.6,0.8,1.0,1.2,1.4])
    ax[0].set_yticklabels([60,80,100,120,140])
    ax[1].set_xlim([-12,37])
    ax[1].set_ylim([0.6,1.4])
    ax[1].set_xticks([-10,-5,0,5,10,15,20,25,30,35])
    ax[1].set_yticks([0.6,0.8,1.0,1.2,1.4])
    ax[1].set_yticklabels([60,80,100,120,140])
    
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
    fig.savefig(os.path.join(output_loc, 'LTD_final.pdf'), dpi=300)
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
    ax[1,0].set_ylim(0,300)
    ax[1,0].set_ylabel('Membrane Resistance (MΩ)')    
    ax[1,1].set_ylim(0.4,1.6)
    ax[1,1].set_ylabel('Normalized Membrane')
    ax[1,2].set_ylim(0,600)
    ax[1,2].set_ylabel('Holding Current (pA)')
    ax[1,3].set_ylim(0.4,1.6)
    ax[1,3].set_ylabel('Normalized Holding')

    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, 'misc_plasticity_info.pdf'))
    return fig, ax


def domain_plot(data):
    cm = 1/2.54
    neg_data = data.loc[(data.domain == 2) & (data.group == 'neg')]
    pos_data = data.loc[(data.domain == 0) & (data.group == 'pos')]
    plot_data = neg_data.append(pos_data)
    pre_data = plot_data.loc[(plot_data.time < 0) & (plot_data.time > -11)]
    post_data = plot_data.loc[(plot_data.time >  0) & (plot_data.time < 36)]
    
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

def extended_data_fig(data, output_loc):
    all_data = data.loc[(data.time > 19) & (data.time < 26)] #Data for general graph, not necessary to have loc
    all_data = all_data.groupby(['group','cell','pre_post']).mean().reset_index()
    plot_data = data
    #plot_data = data.loc[data.confirmed == 'y'] #Make sure we have data that actually has a cell loc
    plot_data = plot_data.loc[(plot_data.time > 19) & (plot_data.time < 26)] #All data from 20-25 minutes ('long term')
    plot_data = plot_data.groupby(['group','cell','pre_post']).mean().reset_index()
    
    cm=1/2.54
    fig, ax = plt.subplots(1,2, dpi=300, figsize=(0.3*(18.3*cm),3),gridspec_kw={'width_ratios': [1, 2]})
    sns.boxplot(data=all_data, x='group', y='norm_EPSC1_amplitude', fliersize=0, ax=ax[0], order=['pos','neg'], linewidth=0.5)
    sns.boxplot(data=plot_data, x='domain', y='norm_EPSC1_amplitude', hue='group', fliersize=0, ax=ax[1], 
                hue_order=['pos','neg'], linewidth=0.5)
    sns.stripplot(data=all_data, x='group', y='norm_EPSC1_amplitude', ax=ax[0], order=['pos','neg'], size=4, hue='group',
                  hue_order=['pos','neg'], jitter=0.05, marker='^', edgecolor='black', color='white', linewidth=0.5)
    sns.stripplot(data=plot_data, x='domain', y='norm_EPSC1_amplitude', ax=ax[1], hue='group', dodge=True, size=4,
                  hue_order=['pos','neg'], jitter=0.05, marker='^', edgecolor='black', color='white', linewidth=0.5)
    
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.set_ylabel('EPSC (%) (t=20-25 min)')
        a.set_ylim(0.0, 1.2)
        a.set_yticks([0.0,0.4,0.8,1.2])
        
    ax[0].legend().set_visible(False)
    ax[1].legend().set_visible(False)
    ax[0].set_xticklabels(['E+','E-'])
    ax[1].set_xticklabels(['AZ','CZ','PZ/NZ'])
    ax[0].set_xlabel(' ')
    ax[1].set_xlabel(' ')

    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "extendedLTDdata.pdf"), dpi = 300)
    return fig, ax

def extended_data_flocculus(data, output_loc):
    all_data = data.loc[(data.time > 19) & (data.time < 26)] #Data for general graph, not necessary to have loc
    all_data = all_data.groupby(['group','cell','pre_post']).mean().reset_index()
    plot_data = data
    plot_data = plot_data.loc[(plot_data.time > 19) & (plot_data.time < 26)] #All data from 20-25 minutes ('long term')
    plot_data = plot_data.groupby(['group','cell','pre_post']).mean().reset_index()
    #plot_data2 = plot_data.loc[plot_data.group =! 'flocculus']
    cm=1/2.54
    fig, ax = plt.subplots(1,2, dpi=300, figsize=(9.15*cm,3))
    sns.boxplot(data=all_data, x='group', y='norm_EPSC1_amplitude', fliersize=0, ax=ax[0], order=['pos','neg','flocculus'])
    sns.boxplot(data=plot_data, x='domain', y='norm_EPSC1_amplitude', hue='group', fliersize=0, ax=ax[1], 
                hue_order=['pos','neg','flocculus'])
    sns.stripplot(data=all_data, x='group', y='norm_EPSC1_amplitude', ax=ax[0], order=['pos','neg','flocculus'], color='black',
                  jitter=0.05, marker='^')
    sns.stripplot(data=plot_data, x='domain', y='norm_EPSC1_amplitude', hue='group', ax=ax[1], dodge=True, fliersize=0,
                  hue_order=['pos','neg','flocculus'], color='black', jitter=0.05, marker='^')
    
    for a in ax.flatten():
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
        a.set_ylabel('EPSC (%)')
        a.set_ylim(0.3, 1.2)
        #a.set_yticklabels([30,40,50,60,70,80,90,100,110,120])
        
    ax[0].legend().set_visible(False)
    ax[1].legend().set_visible(False)
    ax[0].set_xticklabels(['E+','E-','Flcls'])
    ax[1].set_xticklabels(['AZ','CZ','PZ'])
    ax[0].set_xlabel('EAAT4 Identity')
    ax[1].set_xlabel('Region')

    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "extendedLTDdata_wflocculus.png"), dpi = 1000)
    return fig, ax


def export_data(data):
    #for RMAnova
    post_data = data.loc[data.pre_post == 'post']
    post_data = post_data.loc[post_data.time < 35]
    
    a = np.arange(4,30,5)
    for i in a:
        df = post_data.loc[(post_data.time > i) & (post_data.time <i+6)]
        #print(df)
        print(i+6)
        df = df.groupby(['cell','group']).mean().reset_index()
        df.to_excel('data' + str(i)+'.xlsx')
    return 

if __name__ == "__main__":
    base_folder = "//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data"
    data_loc = os.path.join(base_folder, "Data/LTD/Data")
    meta_loc = os.path.join(base_folder, "Metadata/ltdmeta_select.xlsx")
    output_loc = os.path.join(base_folder, "Analysis_tools/output/")
    plt.style.use('seaborn-paper')
    #color = ["#FFFFFF","#808080"] #white/grey
    color = ["#5BE12D","#D24646"] #Original green/red
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

