import json
import os
import glob
import re

for file in glob.glob("./../stage/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
        gene_id = data["gene_id"]
        ##
        if data["genbank"]["Source"] == "Escherichia coli str. K-12 substr. MG1655":
            data["info"]["documentation"]["source"] = "Escherichia coli MG1655"
            data["info"]["gene_metadata"]["collection_id"] = 4
        if data["genbank"]["Source"] == "Bacillus subtilis subsp. subtilis str. 168": 
            data["info"]["documentation"]["source"] = "Bacillus subtilis 168"
        data["info"]["gene_metadata"]["cloning"]["optimize"] = True
        del data["info"]["gene_metadata"]["optimize"]
        ##
        with open("./../stage/{}/{}.json".format(gene_id,gene_id),"w+") as json_file:
            json.dump(data,json_file,indent=2)

