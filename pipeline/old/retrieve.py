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
import codecs


def freegenes_retrieve():
    
    ## ==========
    ## Parameters
    ## ==========
    
    counter = 0
    base_pairs = 0
    write = input("Write files to system? (ENTER to continue, 'WRITE' to write) : ")
    
    template = ff.Json_load("./configuration/template.json")
    config = ff.FreeGenes_configuration()
    SINGLE_SUBMISSION_ID = config["SINGLE_SUBMISSION_ID"]
    BULK_SUBMISSION_ID = config["BULK_SUBMISSION_ID"]
    PREVIOUS_CSV_SINGLE = config["PREVIOUS_CSV_SINGLE"]
    PREVIOUS_CSV_BULK = config["PREVIOUS_CSV_BULK"]
    
    ## ======
    ## Script
    ## ======
    
    # Bulk data
    current_data_bulk = ff.fill_strip(ff.extract_google_form(BULK_SUBMISSION_ID))
    previous_bulk = ff.fill_strip(pd.read_csv(PREVIOUS_CSV_BULK))
    bulk_data = ff.uniq_data(current_data_bulk, previous_bulk)
    for index, row in bulk_data.iterrows():
        collection_id = ff.NextCollection()
        csv_url = row['CSV file']
        csv_data = ff.csvtext_to_pandas(ff.get_wufoo_textfile(csv_url)).fillna('')
        project_description = row['project_description']
        for index_csv, row_csv in csv_data.iterrows():
            # Setup import into object
            gene_id = ff.NextID() # New ID
            zip_file = zipfile.ZipFile(io.BytesIO(requests.get(row['Genbank files']).content)) # Download zip file
            for name in zip_file.namelist():
                if row_csv['Exact genbank name in zip file'] in name:
                    with zip_file.open(name) as myfile:
                        genbank_file="\n".join(str(myfile.read(), 'utf-8').splitlines()) # Format genbank file all nice
            ## IMPLEMENT CONDITIONAlS
            try:
                conditionals = row_csv["conditionals"]
            except:
                conditionals = []
            try:
                hashtags = row_csv['Hashtags'].split("#")
            except:
                hashtags = []
            try:
                optimize = row_csv["optimize"]
            except:
                optimize = False
            # Import into object
            freegene = ff.FreeGene(gene_id, collection_id, row['Timestamp'], row['Name'], row['Email'], row['Affiliation'], row['ORCID'], row_csv['Gene name'], row_csv['Description'], row_csv['Links'], row_csv['Part type'], row_csv['Source organism'], row_csv['Target organism'], row_csv['Safety information'], genbank_file, template, optimize, project_description=project_description, tags=[], hashtags=hashtags, conditionals=conditionals)
            counter += 1
            if write == "WRITE":
                freegene.json_write()
            base_pairs = base_pairs + len(freegene.optimized)
    if write == "WRITE":
        current_data_bulk.to_csv(PREVIOUS_CSV_BULK, index=False)
    
    # Single submissions
    collection_id = ff.NextCollection()
    current_data_single = ff.fill_strip(ff.extract_google_form(SINGLE_SUBMISSION_ID))
    previous_single = ff.fill_strip(pd.read_csv(PREVIOUS_CSV_SINGLE))
    single_data = ff.uniq_data(current_data_single, previous_single)
    for index, row in single_data.iterrows():
        gene_id = ff.NextID() # New ID
        genbank_file = ff.get_wufoo_textfile(row['genbank_file'])
        if row['optimize'] == "Yes":
            optimize = "custom_1"
        else:
            optimize = False
        freegene = ff.FreeGene(gene_id, collection_id, row['Timestamp'], row['name'], row['email'], row['affiliation'], row['orcid'], row['gene_name'], row['description'], row['links'], row['part_type'], row['source_organism'],  row['target_organism'], row['safety_info'], genbank_file, template, optimize, project_description=row["project_description"], tags=row["tags"].split(","), hashtags=row["hashtags"].split("#"), conditionals=row["conditionals"])
        counter += 1
        if write == "WRITE":
            freegene.json_write()
        base_pairs = base_pairs + len(freegene.optimized)
    if write == "WRITE":
        current_data_single.to_csv(PREVIOUS_CSV_SINGLE, index=False)
    
    # Text output
    print("Total of " + str(counter) + " genes to synthesize")
    print("Numbering approximately " + str(base_pairs) + " base pairs")
    return freegenes_manager()
    
    
    
