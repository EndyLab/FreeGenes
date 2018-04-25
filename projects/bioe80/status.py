import json
import os
import glob
import re
import csv


for file in glob.glob("./data/*/*.json"):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    # Gene name
    info = {}
    info["name"] = data["info"]["documentation"]["gene_name"]
    # Gene ID
    for ref in data["genbank"]["references"]:
        gene_id = ''.join(re.findall(r'\"GeneID:([A-Za-z0-9:_/-_\s-]+)\"',ref))
        if gene_id:
            info["ncbi"] = 'https://www.ncbi.nlm.nih.gov/gene?term=' + gene_id
    # E.coli wiki
    info["wiki"] = 'http://ecoliwiki.net/colipedia/index.php/' + data["genbank"]["locus_tag"]
    # Protein database
    info["protein"] = 'https://www.ncbi.nlm.nih.gov/protein/' + data["genbank"]["protein_id"]
    # Enzyme database
    if data["genbank"]["EC_number"]:
        info["enzyme"] = 'https://enzyme.expasy.org/EC/' + data["genbank"]["EC_number"]
    else:
        info["enzyme"] = ""
    print(info)
    with open('spreadsheet.csv','a') as f:
        w = csv.writer(f)
        w.writerow(info.values())

