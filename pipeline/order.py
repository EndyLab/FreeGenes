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

config = ff.FreeGenes_configuration()
stage = config["STAGE_PATH"]

random_sequence = "CGTAACTCGATCACTCACTC"
small_seq_ids = []
small_seqs = []
large_seq_ids = []
large_seqs = []


def fragment_genes():
    for file in glob.glob(stage + "*/*.json"):
        with open(file,"r") as json_file:
            data = json.load(json_file)
            part_type = data["info"]["gene_metadata"]["cloning"]["part_type"]
            part_type = part_type.lower()
            part_type = part_type.replace(" ", "_")
            if not part_type == "vector":
                ortho_pair = ff.get_pair(data["gene_id"])
            else:
                ortho_pair = ("","")
            fragments = ff.FG_standard_fragment(data["sequence"]["optimized_sequence"],part_type,data["info"]["gene_metadata"]["cloning"]["cloning_enzyme"], ortho_pair)
            print(fragments)
            print(data["gene_id"])
            gene_id = data["gene_id"]
            for index,frag in enumerate(fragments):
                fragment_name = gene_id + "_" + str(index + 1)
                data["sequence"]["fragment_sequences"][fragment_name] = frag
            path = "{}{}".format(stage,gene_id)
            with open("{}/{}.json".format(path,gene_id),"w+") as json_file:
                json.dump(data,json_file,indent=2)


def write_link():
    for file in glob.glob(stage + "*/*.json"):
        with open(file,"r") as json_file:
            data = json.load(json_file)
        # Fragment the genes
        # Begin query
        if data["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] == "BtgZI":
            small_seq_ids.append(data["gene_id"])
            small_seqs.append(data["sequence"]["fragment_sequences"]["{}_1".format(data["gene_id"])])
        elif data["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] == "BbsI":
            print("Num frags: ",len(data["sequence"]["fragment_sequences"]))
            if len(data["sequence"]["fragment_sequences"]) > 1:
                print("too many frags")
                continue
            large_seq_ids.append(data["gene_id"])
            large_seqs.append(data["sequence"]["fragment_sequences"]["{}_1".format(data["gene_id"])])
    
    # Generate dataframes that are sorted in opposite directions based on length
    # which pairs the smallest large fragment with the largest small fragment
    small_df = pd.DataFrame({
        "Gene ID" : small_seq_ids,
        "Sequence" : small_seqs,
        "Length" : [len(seq) for seq in small_seqs]
    })
    small_df = small_df.sort_values("Length",ascending=False)
    large_df = pd.DataFrame({
        "Gene ID" : large_seq_ids,
        "Sequence" : large_seqs,
        "Length" : [len(seq) for seq in large_seqs]
    })
    large_df = large_df.sort_values("Length")
    
    small_counter = 0
    print("Total small sequences: ",len(small_df))
    
    ## ====================================================
    ## Join Fragments
    ## ====================================================
    joined_seqs = []
    joined_ids = []
    fragment_names = []
    
    # Pair sequences together until it runs out of either type of sequence
    for index,row in large_df.iterrows():
        print("small counter: ",small_counter)
        if len(small_df) == small_counter:
            print("ran out of small")
            break
        small_row = small_df.iloc[small_counter]
        if ff.repeat_check(row["Sequence"] + small_row["Sequence"]):
            joined_seq = row["Sequence"] + small_row["Sequence"]
            joined_ids.append(row["Gene ID"])
            joined_seqs.append(joined_seq)
            fragment_names.append(row["Gene ID"] + "_link_" + small_row["Gene ID"])
            joined_ids.append(small_row["Gene ID"])
            joined_seqs.append(joined_seq)
            fragment_names.append(row["Gene ID"] + "_link_" + small_row["Gene ID"])
            small_counter += 1
    joined_df = pd.DataFrame({
        "Gene ID" : joined_ids,
        "Sequence" : joined_seqs,
        "Fragment Name" : fragment_names
    })
    
    # Change the files in the database to reflect the joined sequences
    for index,row in joined_df.iterrows():
        with open("{}{}/{}.json".format(stage,row["Gene ID"],row["Gene ID"]),"r") as json_file:
            data = json.load(json_file)
        data["sequence"]["fragment_sequences"] = {}
        data["sequence"]["fragment_sequences"][row["Fragment Name"]] = row["Sequence"]
        with open("{}{}/{}.json".format(stage,row["Gene ID"],row["Gene ID"]),"w+") as json_file:
            json.dump(data,json_file,indent=2)

def twist_order():
    next_sub_num = input("Next submission number : ")
    ## Find all of the sequences that have yet to be ordered
    will_order = []
    will_order_seqs = []
    for file in glob.glob("{}*/*.json".format(stage)):
        with open(file,"r") as json_file:
            data = json.load(json_file)
    
        # Excludes sequences that have already been ordered and small sequences
        # that haven't been paired yet
        #if data["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] == "BtgZI":
            #continue
        # Only pulls the sequence to order from the large fragment
        if data["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] == "BbsI":
            for fragment in data["sequence"]["fragment_sequences"]:
                print("fragment",fragment)
                will_order.append(fragment)
                will_order_seqs.append(data["sequence"]["fragment_sequences"][fragment]) 
        data["info"]["order_number"] = int(next_sub_num)
        with open(file,"w+") as json_file:
            json.dump(data,json_file,indent=2)
    
    # Output DNA in Twist order format
    twist_dna = pd.DataFrame({
            'gene name': will_order,
            'FASTA_seq': will_order_seqs,
            }, columns=['gene name','FASTA_seq']) 
    
    previous_submissions = (sorted(glob.glob("./.." + "/submissions/*.csv")))
    twist_dna.to_csv('{}/submissions/submission{}.csv'.format("./..",str(next_sub_num).zfill(3)),index=False)
    print("Completed submission form.")

def replace_bad_sequence(gene_id):
    table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)
    json_data = ff.Json_load(stage + "{}/{}.json".format(gene_id,gene_id))
    seq_to_replace = input("Sequence to replace? ")
    new_seq = input("New sequence? ")
    old_sequence = json_data["sequence"]["optimized_sequence"]
    new_sequence = ff.fix_sequence(table,gene_id,old_sequence.replace(seq_to_replace, new_seq),0,[])
    if not Seq(old_sequence, IUPAC.unambiguous_dna).translate() == Seq(new_sequence, IUPAC.unambiguous_dna).translate():
        print("Bad translation, try again")
        replace_bad_sequence(gene_id)
    else:
        with open(stage + "{}/{}.gb".format(gene_id,gene_id),"r") as genbank_single:
            genbank_current = genbank_single.read()
        genbank_fixed = ff.replace_genbank_sequence(genbank_current, new_sequence)
        with open(stage + "{}/{}.gb".format(gene_id,gene_id),"w+") as genbank_single:
            genbank_single.write(genbank_fixed)
        json_data["sequence"]["optimized_sequence"] = new_sequence
        with open(stage + "{}/{}.json".format(gene_id,gene_id),"w+") as json_file:
            json.dump(json_data,json_file,indent=2)
        print("Wrote new sequence for " + gene_id)

        if input("Replace another sequence? Y or N : ").upper() == "Y":
            new_gene_id = input("gene id : ")
            return replace_bad_sequence(new_gene_id)
        else: 
            fragment_to_order()

def reoptimize_fragment(gene_id):
    json_data = ff.Json_load(stage + "{}/{}.json".format(gene_id,gene_id))
    table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)
    translation = Seq(json_data["sequence"]["optimized_sequence"], IUPAC.unambiguous_dna).translate()[:-1] 
    #translation = json_data["genbank"]["translation"]
    new_sequence = ff.fix_sequence(table,gene_id,codon.optimize_protein(table, translation) + "TGA")
    if not Seq(json_data["sequence"]["optimized_sequence"], IUPAC.unambiguous_dna).translate() == Seq(new_sequence, IUPAC.unambiguous_dna).translate():
        print("Bad translation, try again")
        reoptimize_fragment(gene_id)
    else:
        with open(stage + "{}/{}.gb".format(gene_id,gene_id),"r") as genbank_single:
            genbank_current = genbank_single.read()
        genbank_fixed = ff.replace_genbank_sequence(genbank_current, new_sequence)
        with open(stage + "{}/{}.gb".format(gene_id,gene_id),"w+") as genbank_single:
            genbank_single.write(genbank_fixed)
        json_data["sequence"]["optimized_sequence"] = new_sequence
        with open(stage + "{}/{}.json".format(gene_id,gene_id),"w+") as json_file:
            json.dump(json_data,json_file,indent=2)
        print("Wrote new sequence for " + gene_id)


def reset_fragment_stage():
    for file in glob.glob(stage + "*/*.json"):
        with open(file,"r") as json_file:
            data = json.load(json_file)
            data["sequence"]["fragment_sequences"] = {}
            gene_id = data["gene_id"]
            path = "{}{}".format(stage,gene_id)
            with open("{}/{}.json".format(path,gene_id),"w+") as json_file:
                json.dump(data,json_file,indent=2)
    print("Fragments on stage cleared.")

def id_reset():
    number_of_files = len(glob.glob(stage + "*")) 
    input("Are you sure you'd like to continue? (ENTER TO CONTINUE)")
    original_number = config["ID_START"]
    new_number = config["ID_START"] + number_of_files
    config["ID_START"] = new_number
    # Collection reset
    config["LAST_COLLECTION"] = config["LAST_COLLECTION"] + 1
    with open("./configuration/FreeGene_config.yaml","w+") as yaml_file:
        yaml.dump(config,yaml_file,default_flow_style=False)
    print("Replaced {} with {}. New collection number start is {}".format(original_number,new_number,config["LAST_COLLECTION"]))


def fragment_to_order():
    print("Recreating database")
    reset_fragment_stage()
    fragment_genes()
    write_link()
    twist_order()


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

