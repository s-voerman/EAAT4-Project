import pandas as pd
import numpy as np

data_folder = "/home/bastiaan/Neuroscience/EAAT4/output/ipsc/"

os.chdir(data_folder)

def grab_data(data_folder):
    output = pd.DataFrame()
    for x in os.listdir(data_folder):
        o = pd.read_excel(x)
        output = output.append(o)
    output = output.reset_index()
    return output

def clean(d):
    data = pd.DataFrame(d)
    data["fileloc"] = data.loc[:, "cell"]
    data["cell"] = [x.split("/")[-1].split(".")[0].split("-")[0] for x in data["fileloc"]]
    data["trace"] = [x.split("/")[-1].split(".")[0].split("-")[1] for x in data["fileloc"]]
    return data

def read_meta():
    meta = pd.read_csv("/home/bastiaan/Neuroscience/EAAT4/data/metadata/zebrin.csv", sep="\t")
    meta.loc[:, "name"] = ["_".join([x.lower(), y]) for x,y in zip(meta.patcher, meta.cell)]
    return meta

def find_location(meta, data):
    for cell in data.cell.unique():
        try:
            zebrin = meta.loc[meta.name == cell, "zebrin"].iloc[0]
            location = meta.loc[meta.name == cell, "lobule"].iloc[0]
        except IndexError:
            zebrin, location = np.nan, np.nan
        data.loc[data.cell == cell, "zebrin"] = zebrin
        data.loc[data.cell == cell, "location"] = location
    data = data.loc[~data.zebrin.isna()]
    return data

def remove_bad_events(data):
    pass

def create_fig(data):
    fig, ax = plt.subplots(4, 1)
    indep = "location"
    counts = data.groupby(["cell", "zebrin", "location"]).count().reset_index()
    medians = data.groupby(["cell", "zebrin", "location"]).median().reset_index()
    sns.barplot(x="location", y="amplitude", data=medians, ci=False, alpha=.5,
                ax=ax[0])
    sns.stripplot(x="location", y="amplitude", data=medians, alpha=1, ax=ax[0])
    sns.barplot(x="location", y="rise", data=medians, ci=False, alpha=.5,
                ax=ax[1])
    sns.stripplot(x="location", y="rise", data=medians, alpha=1, ax=ax[1])
    sns.barplot(x="location", y="tau", data=medians, ci=False, alpha=.5,
                ax=ax[2])
    sns.stripplot(x="location", y="tau", data=medians, alpha=1, ax=ax[2])
    sns.barplot(x="location", y="amplitude", data=counts, ci=False, alpha=.5,
                ax=ax[3])
    sns.stripplot(x="location", y="amplitude", data=counts, alpha=1, ax=ax[3])
    ax[3].set_ylabel("Counts")
    return fig, ax 

def create_kde(data):
    fig, ax = plt.subplots(3, 1)
    zebpos = data.loc[data.zebrin == "pos"]
    zebneg = data.loc[data.zebrin == "neg"]
    bw = [.1, 0.001, 0.001]
    for i, v in enumerate(["amplitude", "tau", "rise"]):
        sns.kdeplot(zebpos.loc[:, v], color="black", ax=ax[i], bw=bw[i], label="z+")
        sns.kdeplot(zebneg.loc[:, v], color="red", ax=ax[i], bw=bw[i], label="z-")
    ax[0].set_xlim(0, 200)
    return fig, ax 

