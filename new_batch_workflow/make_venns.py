import os
from glob import glob
import pandas as pd
from matplotlib_venn import venn2, venn3, venn3_circles
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import random


def to_pretty_name(filename):
    return filename.split("/")[-1]

def build_venn2(sets, a, b, title, ax=None):
    A = len(sets[a])
    B = len(sets[b])
    AB = len(sets[a].intersection(sets[b]))
    A -= AB
    B -= AB

    axx = ax
    if axx is None:
        plt.figure(figsize=(6,4))
        axx = plt.gca()

    # A B AB, C, AC, BC, ABC
    v = venn2(subsets=(1,1,1,1,1,1,1), set_labels=(a,b), ax=axx)
    v.get_label_by_id('10').set_text("\n".join(sorted(sets[a].difference(sets[b]))))
    v.get_label_by_id('01').set_text("\n".join(sorted(sets[b].difference(sets[a]))))
    v.get_label_by_id('11').set_text("\n".join(sorted(sets[a].intersection(sets[b]))))
    axx.set_title(title)

    if ax is None:
        plt.show()
    return AB

# A = 4227
# B = 5575
# AB = 4137
# A -= AB
# B -= AB
# v = venn2(subsets=(1,1,1,1,1,1,1), set_labels=("old","new"))
# v.get_label_by_id('10').set_text(str(A))
# v.get_label_by_id('01').set_text(str(B))
# v.get_label_by_id('11').set_text(str(AB))
# plt.title("Old vs New Species")
# plt.show()
# exit(-1)


os.chdir("..")

print("Hello World")
# tsv_glob = "results/wol/linear_regression/tsv/MS_HHC_*_old.tsv"
# tsv_glob = "results/wol/linear_regression/tsv/MS_HHC_*_new.tsv"
# tsv_glob = "results/wol/linear_regression/tsv/UntreatedMS_HHC_*_old.tsv"
# tsv_glob = "results/wol/linear_regression/tsv/UntreatedMS_HHC_*_new.tsv"
# tsv_glob = "results/wol/linear_regression/tsv/MS_HHC_*.tsv"
# tsv_glob = "results/wol/linear_regression/tsv/MS_HHC_(old|new).tsv"
# tsv_glob = "results/wol/linear_regression/tsv/UntreatedMS_HHC_new.tsv"
# names = list(sorted(glob(tsv_glob)))

tsv_glob = "Old Batch vs New Batch"
names = list(sorted([
    "results/wol/linear_regression/tsv/MS_HHC_old_bonferroni_Bootstrap.tsv",
    "results/wol/linear_regression/tsv/MS_HHC_new_bonferroni_Bootstrap.tsv"
]))
# names = list(sorted([
#     "results/wol/linear_regression/tsv/MS_HHC_old_wald.tsv",
#     "results/wol/linear_regression/tsv/MS_HHC_new_wald.tsv"
# ]))
pretty_names = [to_pretty_name(x) for x in names]

top_ks = {}
k = 50
pthresh = None
all_pvals = []
all_bonferroni_pvals = []
all_dfs = []
for n in names:
    df = pd.read_csv(n, sep="\t")
    df = df.sort_values(by="Pr_diseaseMS")
    corrections = ["bonferroni", "fdr"]
    cols = ["taxonomy", "Coef_diseaseMS", "Pr_diseaseMS"]
    for c in corrections:
        if c + "_diseaseMS" in df.columns:
            cols.append(c + "_diseaseMS")

    all_dfs.append(df[cols])

    if k is not None and pthresh is None:
        top_k = set(df["taxonomy"].iloc[:k].to_list())
    else:
        top_k = set(df[df["Pr_diseaseMS"] < pthresh]["taxonomy"].to_list())
    # random_k = set(random.sample(df["taxonomy"].to_list(), k))
    top_ks[to_pretty_name(n)] = top_k
    all_pvals += df["Pr_diseaseMS"].to_list()

all_dfs = [d.set_index("taxonomy") for d in all_dfs]
combined = all_dfs[0].join(all_dfs[1], lsuffix="new", rsuffix="old")
print(combined)
plt.scatter(combined["Coef_diseaseMSold"], combined["Coef_diseaseMSnew"])
for rid, row in combined.iterrows():
    for s in ["Akkermansia", "Eisenbergiella", "Faecalibacterium", "Ruthenibacterium", "Hungatella"]:
        if rid.startswith(s):
            plt.scatter(row["Coef_diseaseMSold"], row["Coef_diseaseMSnew"], c="red")
            plt.text(row["Coef_diseaseMSold"], row["Coef_diseaseMSnew"], rid)


plt.title("Coef")
plt.xlabel("Coef Disease MS (old)")
plt.ylabel("Coef Disease MS (new)")
plt.axis("equal")
plt.show()

plt.scatter(combined["Pr_diseaseMSold"], combined["Pr_diseaseMSnew"])
for rid, row in combined.iterrows():
    for s in ["Akkermansia", "Eisenbergiella", "Faecalibacterium", "Ruthenibacterium", "Hungatella"]:
        if rid.startswith(s):
            plt.scatter(row["Pr_diseaseMSold"], row["Pr_diseaseMSnew"], c="red")
            plt.text(row["Pr_diseaseMSold"], row["Pr_diseaseMSnew"], rid)

plt.title("Pr")
plt.xlabel("Pr Disease MS (old)")
plt.ylabel("Pr Disease MS (new)")
plt.axis("equal")
plt.show()
for c in corrections:
    if c + "_diseaseMSold" in combined.columns:
        plt.scatter(combined[c+"_diseaseMSold"], combined[c+"_diseaseMSnew"])
        for rid, row in combined.iterrows():
            print(row[c+"_diseaseMSold"], row[c+"_diseaseMSnew"], rid)
            plt.text(row[c+"_diseaseMSold"], row[c+"_diseaseMSnew"], rid)
        plt.title("Significant across batches by " + c)
        plt.xlabel(c + "Disease MS (old)")
        plt.ylabel(c + "Disease MS (new)")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim(right=0.05)
        plt.ylim(top=0.05)
        plt.show()

plt.hist(all_pvals, bins=50)
plt.xlabel("p")
plt.ylabel("count")
plt.title("Histogram of P values")
plt.show()


xm = 1
ym = 1

fig, axes = plt.subplots(xm, ym)
x = 0
y = 0
avg_shared = 0
for i in range(len(pretty_names)):
    if i == 1:
        continue
    if xm == 1 and ym == 1:
        ax = axes
    else:
        ax = axes[y, x]
    x += 1

    print(pretty_names[i])
    print(top_ks[pretty_names[i]])
    print(pretty_names[1])
    print(top_ks[pretty_names[1]])
    shared_count = build_venn2(top_ks, pretty_names[1], pretty_names[i], "", ax=ax)
    avg_shared += shared_count

    if x == xm:
        x = 0
        y += 1
    if y == xm:
        x = 0
        y = 0
        plt.suptitle(tsv_glob + "\nTop " + str(k) + " Consistency")
        plt.show()

print("Average Shared:", avg_shared / 9)

plt.suptitle(tsv_glob + "\nTop " + str(k) + " Consistency")
plt.show()

print("")
