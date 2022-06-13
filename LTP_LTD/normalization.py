import pandas as pd
import numpy as np
import seaborn as sns


#DATA = pd.read_excel("./output/LTP_output.xlsx")

def normalization(data):
    out_pd = pd.DataFrame()
    for cell in data.cell.unique():
        sel = data.loc[data.cell == cell]
        sel = sel.loc[sel.pre_post == "pre"]
        sel.loc[:, "ratio"] = sel.loc[:, "r_cap"] / sel.loc[:, "r_input"]
        for v in ["EPSC1_amplitude", "EPSC2_amplitude", "r_input",
                  "r_membrane", "r_cap", "ratio"]:
            m = sel.loc[:, v].mean()
            sel.loc[:, "norm_"+v] = [z / m for z in sel[v]]
        out_pd = pd.concat([out_pd, sel])
    out_pd = out_pd.loc[out_pd.norm_r_membrane < 5]
    out_pd = out_pd.loc[out_pd.norm_r_membrane > -5]
    return out_pd



