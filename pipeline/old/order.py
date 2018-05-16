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
from zipfile import ZipFile
from io import BytesIO
import requests
import zipfile
import io
import freegenes_functions as ff
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import yaml
import codon

def order_manager():
    print("\n")
    print("=== FreeGenes Order manager ===")
    options = ("Frag -> Order", "Set ID starting point", "Fragment genes in stage", "Write linkers to stage", "Create Twist submission spreadsheet", "Change part of sequence", "Reoptimize single fragment", "Reset fragment stage", "Clear ortho pairs", "Exit")
    choice = ff.option_list(options)
    
    if choice == "Frag -> Order":
        fragment_to_order()
    elif choice == "Set ID starting point":
        id_reset()
    elif choice == "Fragment genes in stage":
        fragment_genes()
    elif choice == "Write linkers to stage":
        write_link()
    elif choice == "Create Twist submission spreadsheet":
        twist_order()
    elif choice == "Change part of sequence":
        replace_bad_sequence(input("gene_id : "))
    elif choice == "Reoptimize single fragment":
        reoptimize_fragment(input("gene_id : "))
    elif choice == "Reset fragment stage":
        reset_fragment_stage()
    elif choice == "Exit":
        sys.exit()
    elif choice == "Clear ortho pairs":
        ff.clear_pairs()

    print("Returning to Order manager")
    return order_manager()

## ============
## Run the code
## ============
print("Welcome to the FreeGenes Order manager")

order_manager()

