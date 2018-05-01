from Bio import SeqIO
from io import StringIO
import requests
import json
import pandas as pd
import re
import os
import glob
import datetime
import shutil 
import sys
import freegenes_functions as ff
from zipfile import ZipFile
from io import BytesIO
import requests
import zipfile
import io


## ==========
## Parameters
## ==========

counter = 0
base_pairs = 0
write = input("Write files to system? (ENTER to continue, 'WRITE' to write) : ")

template = ff.Json_load("./configuration/template.json")
config = ff.FreeGenes_configuration()
SINGLE_SUBMISSION_ID = config["SINGLE_SUBMISSION_ID"]
PREVIOUS_CSV_SINGLE = config["PREVIOUS_CSV_SINGLE"]

collection_id = ff.NextCollection()

## ======
## Script
## ======

# Recreate bulk data as a thing. 
## Bulk data
#current_data_bulk = fill_strip(extract_google_form(BULK_SUBMISSION_ID))
#previous_bulk = fill_strip(pd.read_csv(PREVIOUS_CSV_BULK))
#bulk_data = uniq_data(current_data_bulk, previous_bulk)
#for index, row in bulk_data.iterrows():
#    csv_url = row['CSV file']
#    csv_data = csvtext_to_pandas(get_wufoo_textfile(csv_url)).fillna('')
#    for index_csv, row_csv in csv_data.iterrows():
#        # Setup import into object
#        gene_id = NextID(counter) # New ID
#        zip_file = zipfile.ZipFile(io.BytesIO(requests.get(row['Genbank files']).content)) # Download zip file
#        for name in zip_file.namelist():
#            if row_csv['Exact genbank name in zip file'] in name:
#                with zip_file.open(name) as myfile:
#                    genbank_file="\n".join(str(myfile.read(), 'utf-8').splitlines()) # Format genbank file all nice 
#        # Genomic preparation (to use 'genbank' file section when parsing)
#        if row_csv["genomic"]:
#            genomic = True
#        else:
#            genomic = False
#        if row_csv["optimize"]:
#            optimize = row_csv["optimize_table"]
#        else: 
#            optimize = False
#        # Import into object
#        freegene = FreeGene(gene_id, NextCollection(), row['Timestamp'], row['Name'], row['Email'], row['Affiliation'], row['ORCID'], row_csv['Gene name'], row_csv['Description'], row_csv['Links'], row_csv['Part type'], row_csv['Source organism'], row_csv['Target organism'], row_csv['Safety information'], genbank_file, template, genomic, optimize, False)
#        # self, gene_id, collection_id, timestamp, author_name, author_email, author_affiliation, author_orcid, gene_name, description, database_links, part_type, source_organism, target_organism, safety, genbank_file, template_json, optimize, genbank_dictionary, tags
#        if write == "WRITE":
#            freegene.write()
#        counter = counter + 1
#        base_pairs = base_pairs + len(freegene.optimized)
#if write == "WRITE":
#    current_data_bulk.to_csv(PREVIOUS_CSV_BULK, index=False)

# Single submissions
current_data_single = ff.fill_strip(ff.extract_google_form(SINGLE_SUBMISSION_ID))
previous_single = ff.fill_strip(pd.read_csv(PREVIOUS_CSV_SINGLE))
single_data = ff.uniq_data(current_data_single, previous_single)
for index, row in single_data.iterrows():
    gene_id = ff.NextID() # New ID
    print(row['genbank_file'])
    genbank_file = ff.get_wufoo_textfile(row['genbank_file'])
    print(genbank_file)
    genomic = False
    if row['optimize'] == "Yes":
        optimize = "custom_1"
    else:
        optimize = False
        # Import into object 
    freegene = ff.FreeGene(gene_id, collection_id, row['Timestamp'], row['name'], row['email'], row['affiliation'], row['orcid'], row['gene_name'], row['description'], row['links'], row['part_type'], row['source_organism'],  row['target_organism'], row['safety_info'], genbank_file, template, optimize)
    # self, gene_id, collection_id, timestamp, author_name, author_email, author_affiliation, author_orcid, gene_name, description, database_links, part_type, source_organism, target_organism, safety, genbank_file, template_json, optimize, genbank_dictionary, tags
    if write == "WRITE":
        freegene.write()
    counter = counter + 1
    base_pairs = base_pairs + len(freegene.optimized)
if write == "WRITE":
    current_data_single.to_csv(PREVIOUS_CSV_SINGLE, index=False)

# Text output
print("Total of " + str(counter) + " genes to synthesize")
print("Numbering approximately " + str(base_pairs) + " base pairs")



