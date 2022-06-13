import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
plt.ion()


data = pd.DataFrame()
for x in os.listdir("./output"):
    p = pd.read_excel("./output/"+x)
    data = data.append(p, ignore_index=True)

data = data.loc[data.amplitude < data.rms * -1.645]
data = data.loc[data.amplitude > -500]
data = data.loc[data.tau > 0.001]
data.tau_pams = data.tau_pams / 1000
data.rise_pams = data.rise_pams / 1000

wt = [30718, 130618, 190718, 250718]
data["geno"] = ["wt" if x in wt else "mutant" for x in data["date"]]

fig, ax = plt.subplots(6,4)
for i, x in enumerate(["amplitude", "rise", "tau", "rise_pams", "tau_pams"]):
    p = data.groupby(["date", "cell"]).mean().reset_index()
    sns.barplot(x="date", y=x, data=p, ax=ax[i,0], alpha=.5, ci=False)
    sns.stripplot(x="date", y=x, data=p, ax=ax[i,0], jitter=True, hue="cell")
    for y in data["date"].unique():
        sns.kdeplot(data.loc[data.date == y, x], ax=ax[i, 1], cumulative=True, label=y)
p = data.groupby(["date", "cell"]).count().reset_index()
sns.barplot(x="date", y="rms", data=p, alpha=.5, ci=False, ax=ax[5,0])
sns.stripplot(x="date", y="rms", data=p, jitter=True, ax=ax[5,0])
ax[5,0].set_ylabel("Event Count")

for i, x in enumerate(["amplitude", "rise", "tau", "rise_pams", "tau_pams"]):
    p = data.groupby(["date", "cell", "geno"]).mean().reset_index()
    sns.barplot(x="geno", y=x, data=p, ax=ax[i,2], alpha=.5, ci=False)
    sns.stripplot(x="geno", y=x, data=p, ax=ax[i,2], jitter=True, hue="cell")
    for y in data["geno"].unique():
        sns.kdeplot(data.loc[data.geno == y, x], ax=ax[i, 3], cumulative=True, label=y)
p = data.groupby(["date", "cell", "geno"]).count().reset_index()
sns.barplot(x="geno", y="rms", data=p, alpha=.5, ci=False, ax=ax[5,2])
sns.stripplot(x="geno", y="rms", data=p, jitter=True, ax=ax[5,2])
ax[5,0].set_ylabel("Event Count")

for a in ax.flatten():
    a.legend_.remove()

plt.tight_layout()
