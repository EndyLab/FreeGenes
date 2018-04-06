import re
import pandas as pd
import sys
import os
from io import StringIO
import subprocess
import textwrap

print(os.system('ls genome_sequences/'))
genome = input("Which genome? : ")

# Process Genbank file into dataframe
csv = subprocess.check_output("python2 gb2tab.py -f CDS genome_sequences/" + genome, shell=True)
df = pd.read_csv(StringIO(str(csv, "utf-8")), sep='\t')
df.columns = ['locus_tag', 'sequence', 'exons', 'description']

# Process a string of important values into a dictionary
def dictionary_builder(tag_list, string):
    return {key: value for (key, value) in map(lambda tag: [tag, ''.join(re.findall(r'{}=\"([A-Za-z0-9:_./-_\s-]+)\"'.format(tag),string))], tag_list)}

important_tags = ['locus_tag', 'gene', 'gene_synonyms' 'product', 'note', 'Source', 'GenBank_acc', 'protein_id', 'EC_number', 'translation']


# Build dictionary into genbank compatible format
def dictionary_to_genbank(dictionary):
    value_list = []
    for key, value in dictionary.items():
        value_list.append(str("/" + key + "=" + '"' + value + '"'))
    return value_list

def genbank_multiline(genbank_list):
    multiline= ''
    for item in genbank_list:
        split_list = textwrap.wrap(item, width=58)
        for item in split_list:
            multiline = multiline + "                     " + item + "\n"
    return multiline#.rstrip()


# Digest the table 
for index, row in df.iterrows():
    string = row["description"]
    data = (dictionary_builder(important_tags, string))
    transl_table = '                     /transl_table='+''.join(re.findall(r'transl_table=([0-9]+)',string))
    references = (re.findall(r'(db_xref=\"[A-Za-z0-9:_/-_\s-]+\")',string))
    multiline = (genbank_multiline((dictionary_to_genbank(data))) + transl_table + "\n" + genbank_multiline(references))


