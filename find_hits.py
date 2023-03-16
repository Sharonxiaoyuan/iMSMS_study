import sys
import pandas as pd

f1 = sys.argv[1]
df = pd.read_csv(f1, sep="\t")

df = df[["Coef_diseaseMS", "Pr_diseaseMS", "ID", "taxonomy"]]
df = df.sort_values(by="Pr_diseaseMS", ascending=True)

print(df.head(20))
