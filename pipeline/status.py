import json
import os
import glob
import re

for file in glob.glob("./../stage/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
        gene_id = data["gene_id"]
        if data["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] == "BtgZI":
            print(gene_id)

