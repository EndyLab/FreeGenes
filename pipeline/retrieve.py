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
import optimize as optimize_sequence
from zipfile import ZipFile
from io import BytesIO
import requests
import zipfile
import io
import optimize
from config import *

## ==========
## Parameters
## ==========

counter = 0
base_pairs = 0
for file in glob.glob("template.json"): 
        with open(file,"r") as template_json:
            template = json.load(template_json)
write = input("Write files to system? (ENTER to continue, 'WRITE' to write) : ")


## =========
## Functions
## =========

def NextCollection():
    data = glob.glob("./../data/*/*.json")
    collection_list = [2]
    for json_file in data:
        with open(json_file,"r") as json_data:
            collection_json = json.load(json_data)
            collection_list.append(collection_json["info"]["gene_metadata"]["collection_id"])
    return(int(max([e for e in collection_list if isinstance(e, int)])) + 1)

def NextID(counter=0):
    number = ID_STARTING_POINT + 1 + counter
    string_number = str(number)
    id_number = (string_number.zfill(6))
    full_id = "BBF10K_" + id_number
    return full_id

def csvtext_to_pandas(csvtext, byte_string=False):
    if byte_string == True:
        csv_file= StringIO(str(csvtext, 'utf-8'))
        data = pd.read_csv(csv_file)
        return data
    else: 
        csv_file= StringIO(csvtext)
        data = pd.read_csv(csv_file)
        return data



def extract_google_form(google_id):
    response = requests.get('https://docs.google.com/spreadsheet/ccc?key=' + google_id + '&output=csv')
    assert response.status_code == 200, 'Wrong status code'
    byte_string=(response.content)
    data = csvtext_to_pandas(byte_string, True)
    #csv_file= StringIO(str(byte_string, 'utf-8'))
    #data = pd.read_csv(csv_file)
    return data

def fill_strip(data):
    data = data.fillna('')
    data = strip_df(data)
    return data

#previous = pd.read_csv('previous_submissions.csv')
def uniq_data(current_data, previous):
    current_data = current_data.fillna('')
    current_data = strip_df(current_data)
    previous = previous.fillna('')
    previous = strip_df(previous) 
    united_data = pd.concat([current_data, previous])
    united_data.fillna('')
    united_data_grouped = united_data.groupby(list(united_data.columns))
    uniq_data_idx = [x[0] for x in united_data_grouped.indices.values() if len(x) == 1]
    uniq_data = united_data.iloc[uniq_data_idx]
    return uniq_data

def url_fixer(url):
    link,file_id = re.match(r'(https:\/\/drive.google.com\/)open\?id=([A-Za-z0-9-_]+)',url).groups()
    file_url = link + "uc?id=" + file_id
    return file_url

def get_google_textfile(url, zip_file=False):
    file_url = url_fixer(url)
    if zip_file == True:
        data = requests.get(file_url)
        return data
    else: 
        data = requests.get(file_url)
        data.encoding = 'utf-8'
        data = data.text
        return data

def get_wufoo_textfile(url, zip_file=False):
    if zip_file == True:
        data = requests.get(url)
        return data
    else:         
        data = requests.get(url)
        data.encoding = 'utf-8'
        data = data.text
        return data

def strip_df(df):
    for column in df:
        new_col = []
        for item in df[column]:
            new_col.append(str(item).strip())
        df[column] = new_col
    return df

## =======
## Classes
## =======
class FreeGene:
    """FreeGene class"""
    def __init__(self, gene_id, collection_id, timestamp, author_name, author_email, author_affiliation, author_orcid, gene_name, description, database_links, part_type, source_organism, target_organism, safety, genbank_file, template_json):
        self.gene_id = gene_id
        self.collection_id = collection_id
        self.submission_timestamp = timestamp
        self.author_name = author_name
        self.author_email = author_email
        self.author_affiliation = author_affiliation
        self.author_orcid = author_orcid
        self.gene_name = gene_name
        self.description = description
        self.database_links = database_links
        self.part_type = part_type.lower()
        self.target_organism = target_organism
        self.safety = safety
        self.genbank_file = genbank_file
        self.original_sequence = str(SeqIO.read(StringIO(genbank_file), "genbank").seq)
        if part_type == "CDS":
            self.optimized = optimize_sequence.optimize_gene(self.original_sequence, self.target_organism)
            self.build_ready = True
        else:
            self.optimized = self.original_sequence
            #if self.find_enzyme in ("BbsI", "BtgZI", "BsaI")
        if len(self.optimized) < 300:
            self.cloning_enzyme = "BtgZI"
        else:
            self.cloning_enzyme = "BbsI"
        self.retrieval_enzyme = "BsaI"
    # Functions
    def find_enzyme(self):
        seq = self.optimized
        def reverse_complement(seq):
            return seq.translate(str.maketrans("ATGC","TACG"))[::-1]
        cut_sites = [
            ("BbsI", "GAAGAC"),
            ("BtgZI", "GCGATG"),
            ("BsaI", "GGTCTC")]
        def single_finder(enzyme_cut,sequence):
            if enzyme_cut in sequence:
                return True
            else:
                return False
        for enzyme in cut_sites:
            if not single_finder(enzyme[1],seq) and not single_finder(reverse_complement(enzyme[1]),seq):
                return enzyme[0]
                break
    def write(self): 
        path = "./../data/{}".format(gene_id)
        os.makedirs(path)
        template["gene_id"] = self.gene_id
        template["author"]["name"] = self.author_name
        template["author"]["email"] = self.author_email
        template["author"]["affiliation"] = self.author_affiliation
        template["author"]["orcid"] = self.author_orcid
        template["info"]["documentation"]["gene_name"] = self.gene_name
        template["info"]["documentation"]["description"] = self.description
        template["info"]["documentation"]["database_links"] = self.database_links
        template["info"]["gene_metadata"]["cloning"]["part_type"] = self.part_type
        template["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] = self.cloning_enzyme
        template["info"]["gene_metadata"]["cloning"]["retrieval_enzyme"] = self.retrieval_enzyme
        #template["info"]["gene_metadata"]["cloning"]["target_organism"]["organism_name"] = self.target_organism
        template["info"]["gene_metadata"]["safety"] = self.safety
        template["info"]["gene_metadata"]["collection_id"] = self.collection_id
        template["dates"]["submitted"] = self.submission_timestamp
        # Write taxid, target_organism
        template["sequence"]["original_sequence"] = self.original_sequence
        template["sequence"]["optimized_sequence"] = self.optimized
        with open("{}/{}.json".format(path,gene_id),"w+") as json_file:
            json.dump(template,json_file,indent=2)
        with open("{}/{}.gb".format(path,gene_id),"w+") as genbank_single:
            genbank_single.write(self.genbank_file)


## ======
## Script
## ======

# Bulk data
current_data_bulk = fill_strip(extract_google_form(BULK_SUBMISSION_ID))
previous_bulk = fill_strip(pd.read_csv(PREVIOUS_CSV_BULK))
bulk_data = uniq_data(current_data_bulk, previous_bulk)
for index, row in bulk_data.iterrows():
    csv_url = row['CSV file']
    csv_data = csvtext_to_pandas(get_wufoo_textfile(csv_url)).fillna('')
    for index_csv, row_csv in csv_data.iterrows():
        # Setup import into object
        gene_id = NextID(counter) # New ID
        zip_file = zipfile.ZipFile(io.BytesIO(requests.get(row['Genbank files']).content)) # Download zip file
        for name in zip_file.namelist():
            if row_csv['Exact genbank name in zip file'] in name:
                with zip_file.open(name) as myfile:
                    genbank_file="\n".join(str(myfile.read(), 'utf-8').splitlines()) # Format genbank file all nice 
        # Import into object 
        freegene = FreeGene(gene_id, NextCollection(), row['Timestamp'], row['Name'], row['Email'], row['Affiliation'], row['ORCID'], row_csv['Gene name'], row_csv['Description'], row_csv['Links'], row_csv['Part type'], row_csv['Source organism'], row_csv['Target organism'], row_csv['Safety information'], genbank_file, template)
        if write == "WRITE":
            freegene.write()
        counter = counter + 1
        base_pairs = base_pairs + len(freegene.optimized)
if write == "WRITE":
    current_data_bulk.to_csv(PREVIOUS_CSV_BULK, index=False)

# Single submissions
current_data_single = fill_strip(extract_google_form(SINGLE_SUBMISSION_ID))
previous_single = fill_strip(pd.read_csv(PREVIOUS_CSV_SINGLE))
single_data = uniq_data(current_data_single, previous_single)
for index, row in single_data.iterrows():
    gene_id = NextID(counter) # New ID
    genbank_file = get_wufoo_textfile(row['genbank_file'])
        # Import into object 
    freegene = FreeGene(gene_id, NextCollection(), row['Timestamp'], row['name'], row['email'], row['affiliation'], row['orcid'], row['gene_name'], row['description'], row['links'], row['part_type'], row['source_organism'],  row['target_organism'], row['safety_info'], genbank_file, template)
    if write == "WRITE":
        freegene.write()
    counter = counter + 1
    base_pairs = base_pairs + len(freegene.optimized)
if write == "WRITE":
    current_data_single.to_csv(PREVIOUS_CSV_SINGLE, index=False)

# Text output
print("Total of " + str(counter) + " genes to synthesize")
print("Numbering approximately " + str(base_pairs) + " base pairs")



