import pandas as pd
import sys
import xml_get


df = pd.read_csv("info_igem.csv")
#df[['Sequence', 'Type']] = df['Info'].apply(pd.Series)
thing = int(input("range plz: "))
for uses in range(0,thing):
    df = df[df.TotalReuse != uses]

df["Info"] = df["Info"].apply(lambda x: x.replace("('","").replace(")","").replace(" ","").replace("'","").split(","))
df["Seq"] = df["Info"].apply(lambda x: x[0])
df["Type"] = df["Info"].apply(lambda x: x[1])

df["Length"] = df["Seq"].apply(lambda x: len(x))
print("Total bps: {}".format(df["Length"].sum()))
print("Total genes: {}".format(len(df["Length"])))
print("Avg gene length: {}".format(df["Length"].sum()/len(df["Length"])))


df = df.loc[df["Type"] == 'Coding']


df["Length"] = df["Seq"].apply(lambda x: len(x))
print("Total bps: {}".format(df["Length"].sum()))
print("Total genes: {}".format(len(df["Length"])))
print("Avg gene length: {}".format(df["Length"].sum()/len(df["Length"])))





