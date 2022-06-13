import pandas as pd
import numpy as np
import helper_functions as hf
import matplotlib.pyplot as plt
import seaborn as sns

"""
Grab the baseline data for all different potentiation data, calculate the
kinetics of each. Only use data recorded at room temperature, that seems like
it would be reasonable, no? Anyway, ransack all the meta-reading functions from
the other base functions, concatenate into one big file with all the file
locations etc, and then get yeetin.
"""

def combine_all_meta():
    """
    here we make a big ol list with cell id, file location(s), baseline values,
    genos and preferably SR as well?

    general setup of output file:
    cell_id
    geno
    file_loc
    sr?
    """
    def find_meta_1hz(data_loc="./ltp1hz/"):
        out = pd.DataFrame()
        for x in os.listdir(data_loc):
            if os.path.exists(os.path.join(data_loc, x, ".ignore")):
                continue
            meta_path = os.path.join(data_loc, x, "meta.csv")
            meta = pd.read_csv(meta_path, sep="\t")
            meta.loc[:, "date"] = [str(x).zfill(6) for x in meta.date]
            for i, row in meta.iterrows():
                o = {}
                cell = row.owner[0] + "_" + str(row.date).zfill(6) + " " + str(row.cell)
                file_prefix = str(row.date) + " " + str(row.cell)
                group = row.group
                animal = row.animal
                if isinstance(row.pre, str):
                    pre = row.pre.split(" ")
                else:
                    pre = [str(row.pre)]
                pre = [os.path.join(data_loc, x, "".join(["-".join([file_prefix, str(y)]), ".ABF"])) for y in pre]
                o["cell"] = cell
                o["group"] = group
                o["animal"] = animal
                o["pre"] = pre
                o["dataset"] = "ltp1hz"
                o["sr"] = 50000
                out = out.append(o, ignore_index=True)
        return out

    def find_meta_fsk(data_loc="./fsk/"):
        out = pd.DataFrame()
        for x in os.listdir(data_loc):
            meta_path = os.path.join(data_loc, x, "meta.csv")
            meta = pd.read_csv(meta_path, sep="\t")
            meta.loc[:, "date"] = [str(x).zfill(6) for x in meta.date]
            for i, row in meta.iterrows():
                o = {}
                cell = row.owner[0] + "_" + str(row.date).zfill(6) + " " + str(row.cell)
                file_prefix = str(row.date) + " " + str(row.cell)
                group = row.group
                animal = row.animal
                if isinstance(row.pre, str):
                    pre = row.pre.split(" ")
                else:
                    pre = [str(row.pre)]
                pre = [os.path.join(data_loc, x, "".join(["-".join([file_prefix, str(y)]), ".ABF"])) for y in pre]
                o["cell"] = cell
                o["group"] = group
                o["animal"] = animal
                o["dataset"] = "fsk"
                o["pre"] = pre
                o["sr"] = row.sr
                out = out.append(o, ignore_index=True)
        return out

    def find_meta_ltdrt(data_loc="./ltd_rt/"):
        out = pd.DataFrame()
        for x in os.listdir(data_loc):
            if os.path.exists(os.path.join(data_loc, x, ".ignore")):
                continue
            meta_path = os.path.join(data_loc, x, "meta.csv")
            meta = pd.read_csv(meta_path, sep="\t")
            meta.loc[:, "date"] = [str(x).zfill(6) for x in meta.date]
            for i, row in meta.iterrows():
                o = {}
                cell = row.owner[0] + "_" + str(row.date).zfill(6) + " " + str(row.cell)
                file_prefix = str(row.date) + " " + str(row.cell)
                group = row.group
                animal = row.animal
                if isinstance(row.pre, str):
                    pre = row.pre.split(" ")
                else:
                    pre = [str(row.pre)]
                pre = [os.path.join(data_loc, x, "".join(["-".join([file_prefix, str(y)]), ".ABF"])) for y in pre]
                o["cell"] = cell
                o["group"] = group
                o["animal"] = animal
                o["pre"] = pre
                o["dataset"] = "ltdrt"
                o["sr"] = 50000
                out = out.append(o, ignore_index=True)
        return out

    def find_meta_100hz(data_loc="./ltp100hz/"):
        out = pd.DataFrame()
        all_data = pd.DataFrame()
        meta = pd.read_csv(os.path.join(data_loc, "meta.csv"), sep="\t")
        for i, row in meta.iterrows():
            o = {}
            group = row.group
            cell = row.cell
            file_prefix = row.cell
            if isinstance(row.pre, str):
                pre = row.pre.split(" ")
            else:
                pre = [str(row.pre)]
            pre = [os.path.join(data_loc, "".join(["-".join([file_prefix, str(y)]), ".ABF"])) for y in pre]
            o["cell"] = cell
            o["group"] = group
            o["animal"] = "ltp100hz"
            o["pre"] = pre
            o["dataset"] = "ltp100hz"
            o["sr"] = 10000 
            out = out.append(o, ignore_index=True)
        return out

    o = find_meta_1hz()
    o = pd.concat([o, find_meta_fsk()])
    o = pd.concat([o, find_meta_ltdrt()])
    o = pd.concat([o, find_meta_100hz()])

    outlier = [
        "b_051119 1"
    ]
    o = o.loc[~o.cell.isin(outlier)]

    return o

def calculate_kinetics(meta=None):
    """
    We now have a list of all files. Read the file(s), calculate the location
    of the baseline and peak, and start a fittin.
    """

    ev_start = 0.1
    ev_end = 0.2

    out = pd.DataFrame()

    if not meta:
        meta = combine_all_meta()

    for i, row in meta.iterrows():
        data = np.array([])
        for filename in row.pre:
            d, sr = hf.read_abf(filename)
            sr = int(row.sr)
            data = np.concatenate([data, d])

        EVENT_START = int(ev_start * sr)

        if row.dataset == "fsk":
            EVENT_END = int(0.15*sr)
        elif row.dataset == "ltp100hz":
            EVENT_END = int(0.15*sr)
        else:
            EVENT_END = int(0.2*sr)


        windows = np.arange(0, len(data), sr)
        if row.dataset == "ltp100hz":
            windows = np.arange(0, len(data), int(0.95*sr))

        for win in windows:
            d = data[win:win+sr]
            o = {}
            baseline = np.mean(d[EVENT_START-50:EVENT_START-1])
            peak = np.min(d[EVENT_START:EVENT_END])
            amplitude = peak - baseline

            #fig, ax = plt.subplots()
            #ax.plot(np.linspace(0, len(d)/sr, len(d)), d)

            # calculate kinetics
            start = np.where(d == np.max(d[EVENT_START-10:EVENT_END]))[0][0]
            end = np.where(d == peak)[0][0]
            start = np.where(d[start:end] < baseline)[0][0]+start
            fitting_y = d[start:end]
            fitting_x = np.linspace(0, len(fitting_y)/sr, len(fitting_y))
            fit_vars = hf.rise_fit_func(fitting_x, fitting_y, [d[start], 0.01, 0])[0]
            rise = fit_vars[1]*1000

            #ax.plot(fitting_x+(start/sr), hf.rise_func(fitting_x, fit_vars[0], fit_vars[1], fit_vars[2]))

            start = int(end+5)
            end = EVENT_END - 5
            fitting_y = d[start:end]
            print(fitting_y)
            fitting_x = np.linspace(0, len(fitting_y)/sr, len(fitting_y))
            fit_vars = hf.decay_fit_func(fitting_x, fitting_y, [peak, 0.01, 0])[0]
            decay = fit_vars[1]*1000
            #ax.plot(fitting_x+(start/sr), hf.exp_func(fitting_x, fit_vars[0], fit_vars[1], fit_vars[2]))
            #ax.set_xlim(0, 0.2)
            #input()

            # fill in dict
            o["cell"] = row.cell
            o["group"] = row.group
            o["files"] = row.pre
            o["baseline"] = baseline
            o["amplitude"] = amplitude
            o["rise"] = rise
            o["decay"] = decay
            out = out.append(o, ignore_index=True)
    return out

def cleaning_for_paper(data):
    data = data.loc[data.rise < 100]
    data = data.loc[data.decay < 1000]
    data = data.loc[data.rise > 0]
    data = data.loc[data.decay > 0]
    data.loc[:, "group"] = ["wt" if x == "wt" else "mut" for x in data.group]
    data = data.groupby(["cell", "group"]).mean().reset_index()
    print(data.head())

    def pl(data, var, ax):
        sns.barplot(x="group", y=var, data=data, ci=False, ax=ax, alpha=0.5)
        sns.stripplot(x="group", y=var, data=data, size=4, ax=ax)
        return

    fig, ax = plt.subplots(1,3)
    pl(data, "amplitude",ax[0])
    pl(data, "rise", ax[1])
    pl(data, "decay", ax[2])

    return data

