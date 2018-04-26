import re
import pandas as pd
import sys
import os
from io import StringIO
import subprocess
import textwrap
import glob
import json
import datetime
import time
import freegenes_functions as ff
from Bio import SeqIO

## ===============
## SETUP VARIABLES
## ===============

date = datetime.date.today().strftime("%d") + "-" + datetime.date.today().strftime("%B")[:3].upper() + "-" + datetime.date.today().strftime("%Y")
stage = ff.FreeGenes_configuration()["STAGE_PATH"]
counter = 0 

NextCollection = ff.NextCollection()
unused = set()

prefix_genbank = """LOCUS       {}    {} bp ds-DNA     linear   BCT {}
DEFINITION  {}
ACCESSION   .
VERSION     .
KEYWORDS    .
REFERENCE   1  (bases 1 to {})
SOURCE      synthetic DNA sequence
  ORGANISM  {}
  AUTHORS   {}
  TITLE     Direct Submission
FEATURES             Location/Qualifiers
     CDS             {}"""

important_tags = ['locus_tag', 'gene', 'gene_synonyms', 'product', 'note', 'Source', 'GenBank_acc', 'protein_id', 'EC_number', 'translation']

for file in glob.glob("./../pipeline/template.json"):
        with open(file,"r") as template_json:
            template = json.load(template_json)

## =========
## Functions
## =========
# Process a string of important values into a dictionary
def dictionary_builder(tag_list, string):
    return {key: value for (key, value) in map(lambda tag: [tag, ''.join(re.findall(r'/{}=\"([A-Za-z0-9:_./-_\s-]+)\"'.format(tag),string))], tag_list)}

# Build dictionary into genbank compatible format
def dictionary_to_genbank(dictionary):
    value_list = []
    for key, value in dictionary.items():
        if type(value) == type(""):
            value_list.append(str("/" + key + "=" + '"' + value + '"'))
    return value_list

def genbank_multiline(genbank_list):
    multiline= ''
    for item in genbank_list:
        split_list = textwrap.wrap(item, width=58)
        for item in split_list:
            multiline = multiline + "                     " + item + "\n"
    return multiline#.rstrip()

def fasta_refseq_dictionary(file_name):
    dictionary = {}
    for record in SeqIO.parse(file_name, "fasta"):
        dictionary[record.id[:14]] = str(record.seq)
    return dictionary

def genelist_to_list(file_name):
    with open(file_name,"r") as data_file:
        data = data_file.read()
    essential = list(filter(None, data.split("\n")))
    return essential 

## ==================
## Ask for user input
## ==================
path = ff.FreeGenes_configuration()["PROJECT_PATH"]

print("Please choose a genome file")
genome = ff.file_list(path + "genome/genome_sequences/*.gb")
print("Please choose configuration file")
digest_configuration = ff.file_list(path + "genome/digestion_configuration/*.yaml")
config = ff.load_configuration(digest_configuration)
transl_table = '                     /transl_table='+str(config["transl_table"])


## =====================================
## Check for protein files or gene lists
## =====================================
if config["protein_file"]:
    protein_fasta = ff.file_list(path + "genome/protein_lists/*.fasta")
    protein_dictionary = fasta_refseq_dictionary(protein_fasta)

if config["gene_list"]:
    essential_list = ff.file_list(path + "genome/gene_lists/*.txt")
    essential_list = genelist_to_list(essential_list)

## ========================
## Process Genbank into CSV
## ========================
csv_python_command = "python2 " + path+ "genome/gb2tab.py -f CDS "
csv = subprocess.check_output(csv_python_command + genome, shell=True)
df = pd.read_csv(StringIO(str(csv, "utf-8")), sep='\t')
df.columns = ['locus_tag', 'sequence', 'exons', 'description']

## ================
## Digest the table 
## ================
for index, row in df.iterrows():
    string = row["description"]
    sequence = row["sequence"]
    data = (dictionary_builder(important_tags, string))
    references = (re.findall(r'(db_xref=\"[A-Za-z0-9:_/-_\s-]+\")',string))
    data["references"] = references
    # Genbank stuff
    multiline = (genbank_multiline((dictionary_to_genbank(data))) + transl_table + "\n" + "                     /codon_start=1" + "\n" + genbank_multiline(references))
    if not data["gene"] == "":
        definition = data["gene"]
    else:
        definition = data["locus_tag"]

    # Get gene_id and collection id
    gene_id = ff.NextID(counter)
    
    write = True
    # Essential checker
    if config["gene_list"]:
        if data["gene"] in essential_list or data["locus_tag"] in essential_list:
            write = True
        else:
            write = False
            unused.add(definition)

        # Protein file checker
    if config["protein_file"]:
        ref = data["protein_id"]
        if ref:
            write = True
            data["translation"] = protein_dictionary[data["protein_id"]]
            print("Changed translation on " + gene_id)
        else:
            write = False


    # Setup database links
    links = []

    # Setup Genbank file prefix
    genbank_file = prefix_genbank.format(gene_id, str(len(sequence)).rjust(11), date, definition, len(sequence), data["Source"], config["author"], "1.." + str(len(sequence))) + "\n" + multiline
    genbank_file += ff.suffix_genbank(sequence) 
    # Create object
    if write:
        freegene = ff.FreeGene(gene_id, NextCollection, date, config["author"], config["email"], "BioBricks Foundation", "NA", definition, config["description"], links, "CDS", data["Source"], config["target_organism"], config["safety"], genbank_file, ff.Json_load("./configuration/template.json"), config["optimization_table"], data, config["tags"])
        freegene.json_write()
        counter = counter + 1
    else:
        print('skipping ' + definition)


