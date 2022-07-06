import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir('//storage/v/vcl15/ddata/NEUW/dezeeuw/Stijn Voerman/Paper_data/Data/LTD')
data = pd.read_excel('CS_counts.xlsx')

mean_data = data.groupby(['Cell','Group']).mean().reset_index()

sns.barplot(x='Group',y='Spikelets',data=mean_data)
#sns.histplot(data=data,x='Spikelets', hue='Group')