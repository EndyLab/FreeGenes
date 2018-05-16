import json
import os
import glob
import re

for genbank_file in glob.glob("./../stage/*/*.gb"):
    # Read in the file
    with open(genbank_file, 'r') as file :
        filedata = file.read()
    # Replace the target string
    filedata = filedata.replace('db_xref=', '/db_xref=')
    # Write the file out again
    with open(genbank_file, 'w') as file:
        file.write(filedata)
