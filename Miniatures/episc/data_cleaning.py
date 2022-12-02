import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
plt.ion()


"""
Things to implement:
    if fit_coefficient and fit_tau < 0, then not a real event.
"""
def get_data(folder):
    data = pd.DataFrame()
    folder = "./output/"+folder+"/"
    for x in os.listdir(folder):
        d = pd.read_excel(folder + x)
        data = data.append(d, ignore_index=True)
    return data

if __name__ == "__main__":
    folder = input("folder_name? >>>")
    data = get_data(folder)

