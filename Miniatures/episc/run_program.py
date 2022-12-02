import numpy as np
import os
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime
from scipy.signal import (firwin, lfilter, butter, deconvolve, filtfilt,
                          savgol_filter, sosfilt, detrend, resample,
                          medfilt, wiener)
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from statsmodels.robust import mad
import pandas as pd
from neo.rawio import AxonRawIO
from abf_class_v4 import process_ABF
from auxiliary_functions import (convert_to_excel, convert_to_pandas)
import time

plt.ion()

# IMPLEMENT! 
# implement the removal of decay search if going down significantly..
# make rms even more lenient. Like, lower, and larger segments. That should
# help.
# Create output of detection parameters for each run
# more lenient multiplier post-


outliers = []
def create_complete_plot(output_dict, orig_data, filename, output_folder):
    pdo = pd.DataFrame()
    o = output_dict[filename]
    for x in o:
        pdo = pdo.append(o[x], ignore_index=True)

    filename = filename.split("/")[-1][:-4]
    for x in np.arange(0, len(orig_data), 50000):
        start = x / 50000
        end = (x+50000) / 50000
        z = pdo.loc[pdo.time.between(start, end, inclusive=False)]
        time.sleep(0.01)
        fig, ax = plt.subplots(figsize=(15, 5))
        ax.plot(np.linspace(x/50000, (x+50000)/50000, 50000),
                orig_data[x:x+50000], linewidth=.4)
        ax.set_xlim(start, end)
        ax.set_ylim(ax.get_ylim()[0]-30, ax.get_ylim()[1]+30)

        # ALL CRITERIA!
        # z = z.loc[z.exact_amplitude < z.rms * -2]  # NOTE: it is 1 now. should be 1.645
        z = z.loc[np.abs(z.amplitude) > z.rms*2]
        z = z.loc[z.tau > 0.001]
        z = z.loc[z.rise > 0.0004]
        z = z.loc[z.fit_tau < 10000]
        z = z.loc[z.fit_tau > 0]
        z = z.loc[z.tau > z.rise * 0.9]
        z = z.loc[z.tau_pams < -100]

        if z.size > 0:
            for i, event in z.iterrows():
                if event["amplitude"] < event["rms"] * -2:
                    color = "green"
                else:
                    color = "orange"
                xs = event["time"]
                xp = event["rise"] + xs
                xt = event["tau"] + xp
                yb = event["baseline"]
                yp = event["amplitude"]
                yt = yb + (np.exp(-1)*yp)
                fit_x = np.array((event["fit_x"]) / 50000) + xp
                fit_y = event["fit_y"]
                ax.scatter(xs, yb, color=color, s=20)
                ax.scatter(xp, yb + yp, color=color, s=20)
                ax.scatter(xt, yt, color=color, s=20)
                ax.plot(fit_x, fit_y, linewidth=1, color="red")
        ax.set_title(filename)
        base = "./pics/" + output_folder + "/"
        if not os.path.exists(base+filename):
            os.makedirs(base+filename)
        fig.savefig(base + filename + "/" + filename + "_" +
                    str(int(x/50000)).zfill(3) + ".png")
        plt.close()
    return


def go(output_loc, outliers=outliers, files=None, input_pause=False, plot=True):
    for x in files:
        all_files_output = {}
        print(datetime.datetime.now(), x)

        processing = process_ABF(x, output_loc, plot)
        out = processing.run_all()

        print("we here")
        if out == None:
            continue
        elif plot:
            all_files_output[x] = out[0]
            orig_data = out[1]
            decon_data = out[2]
            detection_data = out[3]
            fig = out[4]
            ax = out[5]
        else:
            all_files_output[x] = out

        # create_complete_plot(all_files_output, orig_data, x, output_loc)

        # save space in pandas, remove fit data.
        for ev in all_files_output[x]:
            all_files_output[x][ev]["fit_x"] = None
            all_files_output[x][ev]["fit_y"] = None

        # for debugging etc. 
        if input_pause:
            a = input("type end to stop >>> ")
            if a == "end":
                if plot:
                    return all_files_output, orig_data, decon_data, detection_data, fig, ax, processing
                else:
                    return all_files_output, orig_data, decon_data

        print("now here")
        # output to an excel file.
        pdo = pd.DataFrame()
        for x in all_files_output:
            o = all_files_output[x]
            for y in o:
                y = o[y]
                y["cell"] = x
                pdo = pdo.append(y, ignore_index=True)
        filename = x.split("/")[-1][:-4]
        write = pd.ExcelWriter(output_loc+filename+".xlsx")
        pdo.to_excel(write)
        write.save()

    try:
        return all_files_output, orig_data, decon_data, fig, ax, processing
    except NameError:
        return all_files_output, orig_data, decon_data

def generate_filenames(meta, data_loc):
    for index, row in meta.iterrows():
        seal = "-".join([" ".join([row["name"], str(row["cell"])]), str(row["seal"])])
        seal = os.path.join(data_loc, seal) + ".ABF"
        ipsc = "-".join([" ".join([row["name"], str(row["cell"])]),
                         str(row["ipsc"])])
        ipsc = os.path.join(data_loc, ipsc) + ".ABF"
        meta.loc[index, "seal"] = seal
        meta.loc[index, "ipsc"] = ipsc
    return meta


        
if __name__ == "__main__":
    base_folder = "C:/Users/Gebruiker/Desktop/EAAT4_IPSC_EPSC"
    data_loc = os.path.join(base_folder, "IPSCs/data")
    meta_loc = os.path.join(base_folder, "IPSCs/metadata/mipsc_metafile.xlsx")
    meta = pd.read_excel(meta_loc)
    meta = generate_filenames(meta, data_loc)
    files = meta["ipsc"]
    #output_loc = input("enter output file folder name >>> ")
    output_loc = os.path.join(base_folder, "output", "ipsc")
    print(os.getcwd())
    out = go(output_loc, files=files, input_pause=False, plot=False)
