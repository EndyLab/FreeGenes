from Bio import SeqIO
import os
import glob
import datetime
import sys
import freegenes_fixer
import io
from config import *
from Bio.SeqFeature import SeqFeature, FeatureLocation


prefix_genbank = """LOCUS       {}                 {} bp ds-DNA     linear   BCT {}
DEFINITION  {}
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      synthetic DNA sequence
  ORGANISM  {}
  AUTHORS   {}
  TITLE     Direct Submission
  JOURNAL   FreeGenes object {}
FEATURES             Location/Qualifiers
     source          {}
                     /organism="{}"
                     /mol_type="genomic DNA"
     CDS             {}"""



# Locus_name is just the ID
# raw number of base pairs
#04-APR-2018
#1..753




csv_data = "ecoli_tab.csv"
for index, row in csv_data.iterrows():

