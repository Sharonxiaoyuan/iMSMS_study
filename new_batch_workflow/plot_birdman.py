import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
from scipy import stats
import numpy as np

matplotlib.use("TkAgg")

woltka_metadata = pd.read_csv("../results/wol/woltka_metadata.tsv", sep="\t", index_col="#genome")
high_abundance_gotus = pd.read_csv("./old_1per10000.tsv", index_col="passing")

df_ascs = {}
df_descs = {}


woltka_metadata = pd.read_csv("../results/wol/woltka_metadata.tsv", sep="\t")
woltka_metadata = woltka_metadata.set_index("#genome")

for birdman, birdname in [
    ("../results/wol/birdman_results_imsms_old.tsv", "old"),
    ("../results/wol/birdman_results_imsms_new.tsv", "new")
]:
    df = pd.read_csv(birdman, sep="\t")
    df = df.set_index("Feature")
    df["MS_mean"] = df["C(disease, Treatment('Control'))[T.MS]_mean"]
    df["MS_std"] = df["C(disease, Treatment('Control'))[T.MS]_std"]
    df["MS_Z"] = df["MS_mean"] / df["MS_std"]
    df["MS_p"] = stats.norm.sf(np.abs(df["MS_Z"])) * 2  # *2 for 2 tailed

    # Optionally filter to remove all the crap
    df_filtered = high_abundance_gotus.join(df)
    df = df_filtered

    df["MS_bonferroni"] = np.minimum(df["MS_p"] * df.shape[0], 1)

    df = df.join(woltka_metadata[["genus", "species"]])

    # Sort By
    # sort_by = "MS_mean"
    sort_by = "MS_p"

    df_asc = df.sort_values(sort_by, ascending=True)

    df_asc["rank"] = range(df_asc.shape[0])

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df_asc[["rank", "MS_mean", "species"]])

    df_ascs[birdname] = df_asc


combined_asc = df_ascs["old"].join(df_ascs["new"], lsuffix="_old", rsuffix="_new")
combined_asc = combined_asc.join(woltka_metadata[["genus", "species"]])

x = "MS_Z_old"
y = "MS_Z_new"
title = "MS_Z_filtered_with_bonferroni"
ax = plt.scatter(combined_asc[x], combined_asc[y], alpha=0.25)
for rid, row in combined_asc.iterrows():
    if row["MS_bonferroni_old"] < 0.05 and row["MS_bonferroni_new"] < 0.05: #or \
       # row["genus"] in ["Akkermansia", "Eisenbergiella", "Faecalibacterium", "Ruthenibacterium", "Hungatella", "Lachnospiraceae", "Fusicatenibacter"]: #, "Lactobacillus", "Candidatus", "Brachyspira", "Subdoligranulum", "Succinatimonas", "Prevotella", "Ruminococcus"]:
        plt.scatter(row[x], row[y], c="red")
        plt.text(row[x], row[y], row["species"])
# plt.xscale("log")
# plt.yscale("log")
plt.xlabel(x)
plt.ylabel(y)
# plt.axhline(0.05, linestyle=":")
# plt.axvline(0.05, linestyle=":")

plt.title(title)
plt.savefig("../figures/" + title + ".png")
plt.show()

