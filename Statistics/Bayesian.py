import pandas as pd
import numpy as np
import os
from scipy.stats import levene, ttest_ind
from scipy.stats import linregress
from scipy.stats import tstd

import pymc as pm

#testing bayesian stats in Python

def bayes_factor(e):
    return