import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import helper_functions as hf
import pandas as pd
plt.ion()
pd.options.display.max_columns = 20
pd.options.display.max_rows = 130

def process_cell(object):
    def __init__(self, pre_data_loc, post_data_loc, output_loc):
        self.pre_data_loc = pre_data_loc
        self.post_data_loc = post_data_loc
        self.output_loc = output_loc

    def run(self):
        pre_data, sr = hf.read_abf(self.pre_data_loc)
        post_data, sr = hf.read_abf(self.post_data_loc)
        return pre_data, post_data, sr

if __name__ == "__main__":
    pre = "./data/300816 1-26.ABF"
    post = "./data/300816 1-28.ABF"
    output_loc = "./output/"
    go = process_cell(pre, post, output_loc)
    pre_data, post_data, sr = go.run()
