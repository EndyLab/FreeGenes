import json
import os
import glob
import re


for file in glob.glob("./../data/*/*.json"):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    part = data["info"]["documentation"]["gene_name"]
    print(part + "," + data["gene_id"])
    text_file = open("JCVI.csv","a")
    text_file.write(part + "," + data["gene_id"] + "\n")
    text_file.close()

