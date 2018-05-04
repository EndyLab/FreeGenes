import pandas as pd
import sys
import xml_get


df = pd.read_csv("info_igem.csv")
#df[['Sequence', 'Type']] = df['Info'].apply(pd.Series)

for uses in range(0,3):
    df = df[df.TotalReuse != uses]

df["Length"] = df["Info"].apply(lambda x: len(x.replace("('","").replace("'","").split(",")[0]))
print(df["Length"].sum())
print(len(df))


