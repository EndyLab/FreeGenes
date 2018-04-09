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

### SETUP VARIABLES
date = datetime.date.today().strftime("%d") + "-" + datetime.date.today().strftime("%B")[:3].upper() + "-" + datetime.date.today().strftime("%Y")
author = "FreeGenes Team"
counter = 0
email = "koeng101@gmail.com"
description = "Test"
target_organism = "E.coli"
safety = ""
optimization_table = "custom_1"

print(os.system('ls ./../genome/genome_sequences/'))
genome = input("Which genome? : ")

# open on up template file
for file in glob.glob("./../pipeline/template.json"):
        with open(file,"r") as template_json:
            template = json.load(template_json)


# Process Genbank file into dataframe
csv = subprocess.check_output("python2 ./../genome/gb2tab.py -f CDS ./../genome/genome_sequences/" + genome, shell=True)
df = pd.read_csv(StringIO(str(csv, "utf-8")), sep='\t')
df.columns = ['locus_tag', 'sequence', 'exons', 'description']

# Process a string of important values into a dictionary
def dictionary_builder(tag_list, string):
    return {key: value for (key, value) in map(lambda tag: [tag, ''.join(re.findall(r'{}=\"([A-Za-z0-9:_./-_\s-]+)\"'.format(tag),string))], tag_list)}

important_tags = ['locus_tag', 'gene', 'gene_synonyms', 'product', 'note', 'Source', 'GenBank_acc', 'protein_id', 'EC_number', 'translation']


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

# Genbank prefix and generator for genbank suffix

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





# Digest the table 
for index, row in df.iterrows():
    string = row["description"]
    sequence = row["sequence"]
    data = (dictionary_builder(important_tags, string))
    references = (re.findall(r'(db_xref=\"[A-Za-z0-9:_/-_\s-]+\")',string))
    data["references"] = references
    # Genbank stuff
    transl_table = '                     /transl_table='+''.join(re.findall(r'transl_table=([0-9]+)',string))
    multiline = (genbank_multiline((dictionary_to_genbank(data))) + transl_table + "\n" + "                     /codon_start=1" + "\n" + genbank_multiline(references))
    if not data["gene"] == "":
        definition = data["gene"]
    else:
        definition = data["locus_tag"]
    # Get gene_id and collection id
    gene_id = ff.NextID(counter)

    # Setup database links
    links = []

    # Setup Genbank file prefix
    genbank_file = prefix_genbank.format(gene_id, str(len(sequence)).rjust(11), date, definition, len(sequence), data["Source"], author, "1.." + str(len(sequence))) + "\n" + multiline
    genbank_file += ff.suffix_genbank(sequence)
    #print(genbank_file)
    # Create object
    freegene = ff.FreeGene(gene_id, ff.NextCollection(), datetime.datetime.now(), author, email, "BioBricks Foundation", "NA", definition, description, links, "CDS", data["Source"], target_organism, safety, genbank_file, ff.JsonTemplate(), optimization_table, data)
    counter = counter + 1



# gene_id, collection_id, timestamp, author_name, author_email, author_affiliation, author_orcid, gene_name, description, database_links, part_type, source_organism, target_organism, safety, genbank_file, template_json, genomic, optimize, genbank_dictionary
