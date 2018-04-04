import json
import os
import glob
import re

part_set=set()

for file in glob.glob("./../data/*/*.json"):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    part = data["info"]["gene_metadata"]["cloning"]["part_type"]
    part = part.replace(" ", "_")
    part_set.add(part)
    print(part)
    print(data["gene_id"])

print(part_set)
