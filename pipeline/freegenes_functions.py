import collections
import pickle
import getpass
import smtplib
import datetime
import time
from Bio.SeqUtils import MeltingTemp as mt
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
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from difflib import ndiff
from synbiolib import codon
from Bio.Alphabet import IUPAC
import yaml
from numpy.random import choice
import subprocess

## TODO
# -Implement restriction enzyme chooser
# -Get better at choosing enzymes and such.

# -Make bulk uploading more CLEAR.

# Genomes #
# pGKL1 + pGKL2
# phi29
# lambda
# MS2
# qbeta
# RK plasmid
# F plasmid



## =============
## Configuration
## =============
def load_configuration(config):
    with open(config) as configuration:
        config = yaml.load(configuration)
        return config

def FreeGenes_configuration():
    return load_configuration("./configuration/FreeGene_config.yaml")
config = FreeGenes_configuration()

pairs_file = ""
primer_pairs = pd.DataFrame()

## =============
## DNA sequences
## =============
end_codons = {
    'M': 'ATG',
    'W': 'TGG',
    'F': 'TTT',
    'L': 'CTG',
    'I': 'ATT',
    'V': 'GTG',
    'S': 'TCC',
    'P': 'CCA',
    'T': 'ACC',
    'A': 'GCC',
    'Y': 'TAC',
    'H': 'CAT',
    'Q': 'CAG',
    'N': 'AAC',
    'K': 'AAG',
    'D': 'GAT',
    'E': 'GAG',
    'C': 'TGC',
    'R': 'CGC',
    'G': 'GGC',
    '*': 'TGA'
}

cut_sites= [
        ("BbsI", "GAAGAC"),
        ("BtgZI", "GCGATG"),
        ("BsaI", "GGTCTC"),
        ("BsmBI", "CGTCTC"),
        ("AarI", "CACCTGC"),
        ("BfuAI", "ACCTGC"),
        ("SapI", "GAAGAGC")]

universal_set = pd.read_csv(config["UNIVERSAL_PRIMER"])["Sequence"].tolist()

important_tags = ['locus_tag', 'gene', 'gene_synonyms', 'product', 'note', 'Source', 'GenBank_acc', 'protein_id', 'EC_number', 'translation']

default_table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)


## ===============
## Email functions
## ===============

#gmail_password = getpass.getpass("Enter Email Password: ")
def freegenes_sendmail(gmail_password, send_address, subject, body):
    gmail_user = 'bbfreegenes@gmail.com'

    message = 'From: {}\nTo: {}\nSubject: {}\n\n{}'.format(gmail_user,", ".join(send_address),subject,body)

    try:
        server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
        server.ehlo()
        server.login(gmail_user, gmail_password)
        server.sendmail(gmail_user, send_address, message)
        server.close()

        print ('Email sent!')
    except:
        print ('Something went wrong...')


## ==============
## File Functions
## ==============

def option_list(options):
    return options[int(input(''.join([str(counter) + ". " + str(option) + "\n" for counter, option in list(enumerate(options, start=1))]))) - 1]

def file_list(path):
    counter = 0
    path_number = []
    for files in glob.glob(path):
        counter += 1
        file_name = files.split("/")[-1]
        print("{}. {}".format(counter,file_name))
        path_number.append(files)
    number = input("Which file: ")
    number = int(number) - 1
    return path_number[number]

def Json_load(path):
    for file in glob.glob(path):
        with open(file,"r") as template_json:
            template = json.load(template_json)
        return template

def NextCollection():
    data = glob.glob("./../stage/*/*.json")
    collection_number = config["LAST_COLLECTION"]
    collection_list = [collection_number]
    if data:
        for json_file in data:
            with open(json_file,"r") as json_data:
                collection_json = json.load(json_data)
                collection_list.append(collection_json["info"]["gene_metadata"]["collection_id"])
    return(int(max([e for e in collection_list if isinstance(e, int)])) + 1)

def NextID(count=0):
    counter = len(glob.glob("./../stage/" + "*")) + count
    number = config["ID_START"] + 1 + counter
    string_number = str(number)
    id_number = (string_number.zfill(6))
    full_id = "BBF10K_" + id_number
    return full_id


## ==================
## Retrieve functions
## ==================
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

## ===============
## Digest functions
## ================

def dictionary_builder(tag_list, string):
    return {key: value for (key, value) in map(lambda tag: [tag, ''.join(re.findall(r'/{}=\"([A-Za-z0-9:_./-_\s-]+)\"'.format(tag),string))], tag_list)}

def genbank_to_csv(genome):
    csv_python_command = "python2 " + config["PROJECT_PATH"] + "genome/gb2tab.py -f CDS "
    csv = subprocess.check_output(csv_python_command + genome, shell=True)
    df = pd.read_csv(StringIO(str(csv, "utf-8")), sep='\t')
    df.columns = ['locus_tag', 'sequence', 'exons', 'description']
    return df

# Functions for the digester
def dictionary_to_json_genbank(template, dictionary):
     for key, value in dictionary.items():
         template["genbank"][key] = data[key]
     return template

def suffix_genbank(sequence):
    multiline = "ORIGIN\n"
    start_site = 1
    sequence_list = [sequence[i:i+10] for i in range(0, len(sequence), 10)]
    sequence_split_list = [sequence_list[i:i+6] for i in range(0, len(sequence_list), 6)]
    for index in sequence_split_list:
        multiline += str(start_site).rjust(9)
        for sequence in index:
            multiline += " " + sequence
        multiline += "\n"
        start_site = start_site + 60
    multiline += '//'
    return multiline



## ==================
## Orthogonal primers
## ==================

def load_pairs(file = config["ORTHOGONAL_PRIMER"]):
    global primer_pairs
    global pairs_file
    pairs_file = file
    primer_pairs = pd.read_csv(pairs_file, index_col=0)

def save_pairs():
    primer_pairs.to_csv(pairs_file)

def clear_pairs():
    primer_pairs['GeneID'] = np.nan
    save_pairs()

def get_pair(gene_id, counter=0):
    index = primer_pairs['GeneID'].isnull().idxmax() + counter
    primer_pairs.loc[index, 'GeneID'] = gene_id
    save_pairs()
    return (
        primer_pairs.loc[index, 'Forward'],
        primer_pairs.loc[index, 'Reverse']
    )


### Run pairs
load_pairs()



## =============
## DNA functions
## =============

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATGC","TACG"))[::-1]

# No stop codon on protein or DNA
def transl_checker(DNA, protein):
    if Seq(DNA, IUPAC.unambiguous_dna).translate() == protein:
        return True
    else:
        return False

def FG_MoClo_ends(seq):
    return seq[:-6] + ''.join(list(map(lambda x: end_codons[str(Seq(x, IUPAC.unambiguous_dna).translate())], [seq[-6:-3],seq[-3:]])))

def random_dna_sequence(length):
    return ''.join(np.random.choice(('A', 'C', 'T', 'G')) for _ in range(length))

def translate(seq):
    return str(Seq(seq, generic_dna).translate(table=11))

def sequence_search(search,sequence):
    if search in sequence or reverse_complement(search) in sequence:
        return True
    else:
        return False

def which_vector_enzyme(seq, cut_sites=cut_sites):
    seq = str(seq)
    for enzyme in cut_sites:
        if not sequence_search(enzyme[1],seq):
            return enzyme[0]

def find_enzyme(seq, cut_sites=cut_sites):
    seq = str(seq)
    for enzyme in cut_sites:
        if sequence_search(enzyme[1],seq):
            return enzyme[0]

def max_homopolymer(seq):
    prev = ''
    count = 0
    max_count = 0
    loc = None
    for i, char in enumerate(seq):
        if char == prev:
            count += 1
        else:
            count = 0
        if count > max_count:
            max_count = count
            loc = i - count
        prev = char
    return max_count + 1, loc

## ===============
## Recode_sequence
## ===============

def recode_sequence(table, seq, rep, change_list=[], halfway=True):
    if seq.find(rep) < 0:
        rep = reverse_complement(rep)
        if seq.find(rep) < 0:
            return seq
    # This block contains code that picks a new codon.
    def fraction_list(number_list):
        return list(map(lambda x: x / sum(number_list), number_list))
    def new_codon_chooser(table,old_codons):
        del_table = table
        for codons in old_codons:
            recode_aa = codon_to_aa(table,codons)
            del_table = deletion_table(del_table,recode_aa,codons)
        return ''.join(choice(list(del_table.loc[recode_aa].index), 1, p=fraction_list(list(del_table.loc[recode_aa]['Fraction']))))
    def codon_to_aa(table,old_codon):
        return table.loc[(slice(None),old_codon),:].index[0][0]
    def deletion_table(table,recode_aa,old_codon):
        return table.drop((recode_aa.upper(),old_codon.upper()))

    # This block contains code that produces the desired position information
    def choose_position(list_of_positions, halfway):
        if halfway:
            position = list_of_positions[int(len(list_of_positions)/2)  - 1] 
        else:
            position = list_of_positions[0]
        return position
    def make_position_list(seq, rep):
        position = seq.find(rep)
        position -= position % 3
        position +=1
        list_of_positions = list(range(position, position + (len(rep) // 3) * 3, 3))
        return list_of_positions
    def verify_position(table,old_codon_list):
        new_table = table
        for old_codon in old_codon_list:
            codon_aa = codon_to_aa(table,old_codon)
            new_table = deletion_table(new_table,codon_aa,old_codon)
        to_return = new_table.loc[codon_aa].shape[0]
        return to_return
    def remove_from_list(pos_list,change_list):
        for change in change_list:
            if change in pos_list:
                pos_list.remove(change)
        return pos_list
    def remove_single_aa(table, position_list, change_list):
        # Grab grab spot to modify and its current codon
        position_choice = choose_position(position_list, halfway)
        current_codon = seq[position_choice-1:position_choice+2]
        # Grab all other codons used in that spot
        old_codons = list(map(lambda x: x[0], list(filter(lambda x: x[1] == position_choice, change_list))))
        if not current_codon in old_codons:
            old_codons.append(current_codon)
        if verify_position(table,old_codons) == 0:
            position_list.remove(position_choice)
            return remove_single_aa(table, position_list, change_list)
        else:
            return [position_choice, new_codon_chooser(table,old_codons), current_codon]

    new_recoding = remove_single_aa(table,make_position_list(seq, rep),change_list)
    position_choice = new_recoding[0]
    new_codon = new_recoding[1]
    old_codon = new_recoding[2]
    print("{} -> {} at position {}-{}.".format(old_codon,new_codon,position_choice,position_choice+2))
    return [seq[:position_choice-1] + new_codon + seq[position_choice+2:], position_choice, new_codon]


def recode_region(table,gene_id,seq,rep):
    position = seq.find(rep)
    position -= position % 3
    length_position = len(rep) - (len(rep) % 3)
    recode_sequence = seq[position:position+length_position]
    recode_sequence = codon.optimize_protein(table,translate(recode_sequence))
    final_seq = seq[:position] + recode_sequence + seq[position+length_position:]
    return final_seq



## ======
## Checks
## ======

def repeat_finder(string):
    # Shamelessly stolen from https://codereview.stackexchange.com/questions/63329/finding-the-largest-repeating-substring . Credit to Arvind Padmanabhan.
    l = list(string)
    d = collections.deque(string[1:])
    match = []
    longest_match = []
    while d:
        for i, item in enumerate(d):
            if l[i]==item:
                match.append(item)
            else:
                if len(longest_match) < len(match):
                    longest_match = match
                match = []
        d.popleft()
    return ''.join(longest_match)

def repeat_check(seq):
    repeat_sequence = repeat_finder(seq)
    # Twist's melting temperature algorithm is different than ours
    if len(repeat_sequence) < 20 and int(mt.Tm_Wallace(Seq(repeat_sequence))) < 60:
        return True
    else:
        return False


def cutsite_check(seq, cut_sites=cut_sites):
    seq = str(seq)
    for enzyme in cut_sites:
        if sequence_search(enzyme[1],seq):
            return False
    return True

def cutsite_list(part_type, seq, cut_sites=cut_sites):
    seq = str(seq)
    cut_list = []
    for enzyme in cut_sites:
        if sequence_search(enzyme[1],seq):
            cut_list.append(enzyme[0])
    if not part_type == "vector":
        if "BtgZI" in cut_list and "BbsI" in cut_list or "BsaI" in cut_list:
            raise ValueError("Gene contains {}".format(cut_list))
    return cut_list


def homopolymer_checker(seq):
    homopolymer_max, homopolymer_site = max_homopolymer(seq)
    if homopolymer_max > 6:
        return False
    else:
        return True

def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)
def gc_in_range(seq):
    return (gc_content(seq) > 0.3) & (gc_content(seq) <= 0.65)


def gc_range(seq):
    highest_gc = (0, False)
    lowest_gc = (1, False)
    for start_base in range(len(seq)-50):
        gc_seq = seq[start_base:start_base+50]
        gc_percent = gc_content(gc_seq)
        if gc_percent > highest_gc[0]:
            highest_gc = (gc_percent, gc_seq)
        if gc_percent < lowest_gc[0]:
            lowest_gc = (gc_percent, gc_seq)
    return [highest_gc, lowest_gc]

def gc_50bp_acceptable_range(seq):
    high_low = gc_range(seq)
    average = gc_content(seq)
    return not high_low[0][0] - high_low[1][0] > .51

def gc_range_sequence(seq):
    high_low = gc_range(seq)
    average = gc_content(seq)
    if abs(high_low[0][0] - average) > abs(high_low[1][0] - average):
        return high_low[0][1]
    else:
        return high_low[1][1]


def univerisal_primer_checker(seq):
    for primer in universal_set:
        if primer in seq:
            return False
    return True

def universal_primer_finder(seq):
    for primer in universal_set:
        if primer in seq:
            return primer


## =====
## Fixer
## ===== 

def fix_sequence(table,gene_id,seq,fix_attempts=0,change_list=[], bias="NA"):
    fix_attempts = fix_attempts + 1
    max_attempts = 500
    seq = FG_MoClo_ends(seq)
    if fix_attempts > max_attempts:
        print("{} could not be fixed. Please manually check this sequence.".format(gene_id))
        print("Terminating program")
        sys.exit()
    elif not cutsite_check(seq):
        print("Fixing {} cutsites in {}".format(find_enzyme(seq),gene_id))
        recode = recode_sequence(table, seq, dict((x, y) for x, y in cut_sites)[find_enzyme(seq)], change_list, halfway=True)
        new_seq = recode[0]
        change_list.append([recode[2],recode[1]])
        return fix_sequence(table, gene_id, new_seq, fix_attempts,change_list)
    elif not homopolymer_checker(seq):
        print("Fixing homopolymer stretches in {}".format(gene_id))
        homopolymer_max, homopolymer_site = max_homopolymer(seq)
        homopolymer = seq[homopolymer_site:homopolymer_site+homopolymer_max]
        recode = recode_sequence(table, seq, homopolymer, change_list, halfway=True)
        new_seq = recode[0]
        change_list.append([recode[2],recode[1]])
        return fix_sequence(table, gene_id, new_seq, fix_attempts, change_list)
    elif not repeat_check(seq):
        print("Fixing repeat regions in {}".format(gene_id))
        recode = recode_sequence(table, seq, str(repeat_finder(seq)), change_list, halfway=True)
        new_seq = recode[0]
        change_list.append([recode[2],recode[1]])
        return fix_sequence(table, gene_id, new_seq, fix_attempts, change_list)
    #elif not gc_50bp_acceptable_range(seq):
    #    print("Fixing 50bp of outlier GC in {}".format(gene_id))
    #    gc_outlier = gc_range_sequence(seq)
    #    print(gc_outlier)
    #    new_seq = recode_region(table,gene_id,seq,gc_outlier)
    #    return fix_sequence(table, gene_id, new_seq , fix_attempts, change_list)
    elif not gc_in_range(seq):
        print("Fixing GC content in {}".format(gene_id))
        return fix_sequence(table, gene_id, codon.optimize_protein(table, translate(seq)))
    elif not univerisal_primer_checker(seq):
        print("Fixing internal FG primer regions in {}".format(gene_id))
        return fix_sequence(table, gene_id, universal_primer_finder(seq), fix_attempts, change_list)
    else:
        return check_sequence(gene_id,seq)


## =======
## Checker
## =======

def check_sequence(gene_id,seq,part_type="cds"):
    def is_seq(seq):
        return seq.replace(r'\s +', '').strip() != ""
    def is_triplet(seq):
        return len(seq) % 3 == 0
    def has_start(seq):
        return seq[:3] == "ATG" or "GTG" or "TTG"
    def has_stop(seq):
        return translate(seq)[-1] == "*"
    def has_no_internal_stops(seq):
        return not "*" in translate(seq)[:-1]
    def MoClo_ends(seq):
        return end_codons[translate(seq[-6:-3])] == seq[-6:-3] and seq[-3:] == "TGA"
    def size_in_range(seq):
        return (len(seq) >= 300) & (len(seq) <= 1800)
    def gc_content(seq):
        return (seq.count("G") + seq.count("C")) / len(seq)
    def gc_in_range(seq):
        return (gc_content(seq) > 0.3) & (gc_content(seq) <= 0.65)
    if part_type.lower() == "cds":
        sequence_checks = [is_seq, is_triplet, has_start, has_no_internal_stops, MoClo_ends, homopolymer_checker, gc_in_range, cutsite_check, repeat_check]
    else: 
        sequence_checks = [homopolymer_checker,repeat_check,is_seq, gc_in_range]
    for test in sequence_checks:
        if not test(seq):
            raise ValueError("Gene has failed on {}".format(test.__name__))
    print("{} CHECKED.".format(gene_id))
    return seq

def retrieve_checker(seq):
    check_status = "Clear"
    try:
        check_sequence("gene",seq,"other")
    except Exception as e:
        check_status = e
    return check_status

## ===================
## Refactored fragment
## ===================

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATGCN","TACGN"))[::-1]

def sequence_search(search,sequence):
    if search in sequence or reverse_complement(search) in sequence:
        return True
    else:
        return False

enzymes = {
        "AarI" : {
            "seq" : "CACCTGC",
            "jump" : 4,
            "overhang" : 4
            },
        "BbsI" : {
            "seq" : "GAAGAC",
            "jump" : 2,
            "overhang" : 4
            },
        "BfuAI" : { 
            "seq" : "ACCTGC",
            "jump" : 4,
            "overhang" : 4
            },
        "BsaI" : {
            "seq" : "GGTCTC",
            "jump" : 1,
            "overhang" : 4
            },
        "BsmBI" : {
            "seq" : "CGTCTC",
            "jump" : 1,
            "overhang" : 4
            },
        "BtgZI" : { 
            "seq" : "GCGATG",
            "jump" : 10,
            "overhang" : 4
            },
        "SapI" : { 
            "seq" : "GCTCTTC",
            "jump" : 1,
            "overhang" : 3
            }
        }

FG_part_types = {
        "cds" : {
            "prefix" : "GGTCTCNA",
            "suffix" : "AGAGCTTNGAGACC"
            },
        "eukaryotic_promoter" : {
            "prefix" : "GGTCTCNGGAG",
            "suffix" : "AATGNGAGACC"
            },
        "prokaryotic_promoter" : {
            "prefix" : "GGTCTCNGGAG",
            "suffix" : "TACTNGAGACC"
            },
        "rbs" : { 
            "prefix" : "GGTCTCNTACT",
            "suffix" : "AATGNGAGACC"
            },
        "terminator" : {
            "prefix" : "GGTCTCNGCTT",
            "suffix" : "CGCTNGAGACC"
            },
        "operon" : { 
            "prefix" : "GGTCTCNGGAG",
            "suffix" : "CGCTNGAGACC"
            },
        "cds_aari" : {
            "prefix" : "CACCTGCNNNNGGAG",
            "suffix" : "CGCTNNNNGCAGGTG"
            },
        "cds_operon" : {
            "prefix" : "GGTCTCNA",
            "suffix" : "AGAGCTTNGAGACC"
            },
        "rp_selection" : {
            "prefix": "CACCTGCNNNNGGAG",
            "suffix": "CGCTNNNNGCAGGTG"
            }
        }

def part_type_preparer(part_type, seq, suffix="", prefix=""):
    if part_type == "vector":
        seq = seq + seq[-4:]
        return seq
    if part_type == "custom":
        seq = prefix + seq + suffix
        return seq
    else:
        part = FG_part_types[part_type]
        if len(seq) < 300:
            N_replace = "G"
        else:
            N_replace = "A"
        seq = part["prefix"].replace("N", N_replace) + seq + part["suffix"].replace("N", N_replace)
        return seq
        
# seq is sequence with retrieval prefix, retrieval suffix, and standard flanks already added.
def fragmenter(seq, cloning_enzyme_prefix, cloning_enzyme_suffix, synthesis_max=1500):
    # Setup
    num_frags = len(seq) // synthesis_max + 1
    frag_len = len(seq) // num_frags
    frags = []
    for fragment in range(num_frags):
        frag = seq[max(0, fragment * frag_len -2):min((fragment+1) * frag_len + 2,len(seq))]
        frag = cloning_enzyme_prefix + frag + cloning_enzyme_suffix
        frags.append(frag)
    return frags

def FG_standard_fragment(seq, part_type, cloning_enzyme, ortho_pair):
    random_forward = "CATGCTTGCA"
    random_reverse = "GCTCTGAATA"
    cloning_sites = enzymes[cloning_enzyme]["seq"] + ("N" * enzymes[cloning_enzyme]["jump"])
    cloning_prefix = cloning_sites.replace("N" * cloning_sites.count("N"), random_forward[:cloning_sites.count("N")])
    cloning_suffix = reverse_complement(cloning_sites).replace("N" * cloning_sites.count("N"), random_reverse[:cloning_sites.count("N")])
    return fragmenter("GGAG" + ortho_pair[0] + part_type_preparer(part_type, seq) + reverse_complement(ortho_pair[1]) + "CGCT", cloning_prefix, cloning_suffix)


## ================================
## Twist fix gene by reoptimization
## ================================
def replace_genbank_sequence(original_genbank, optimized_sequence):
    return original_genbank.replace(''.join(re.findall(r'(ORIGIN[A-Za-z0-9:_./-_\s-]+//)', original_genbank)), suffix_genbank(optimized_sequence))

## =======
## Classes
## =======
class FreeGene:
    """FreeGene class"""
    def __init__(self, gene_id, collection_id, timestamp, author_name, author_email, author_affiliation, author_orcid, gene_name, description, database_links, part_type, source_organism, target_organism, safety, genbank_file, template_json, optimize, project_description="", genbank_dictionary={}, tags=[], hashtags=[], conditionals={}):
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
        self.template_json = template_json
        self.project_description = project_description
        self.genbank_dictionary = genbank_dictionary
        self.tags = tags
        self.hashtags = hashtags
        self.conditionals = conditionals
        self.original_genbank = genbank_file
        self.original_sequence = str(SeqIO.read(StringIO(genbank_file), "genbank").seq)
        # Sequence modifications
        print("Running {}".format(self.gene_name))
        if self.part_type == "cds":
            ####### FIX: ONLY DOES CODON CHANGING ACCORDING TO THE E COLI TABLE!!
            if optimize == True:
                print("\n\nOptimizing {}. ( {} )".format(self.gene_id, self.gene_name))
                optimization_table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)
                optimized_sequence = codon.optimize_protein(optimization_table, self.genbank_dictionary["translation"]) + "TAA"
            else:
                print("\n\nNot optimizing {}. ( {} )".format(self.gene_id, self.gene_name))
                optimized_sequence = self.original_sequence
                optimization_table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)

            self.optimized = fix_sequence(optimization_table,self.gene_id,optimized_sequence,0,[])
            print("New seq is {}".format(self.optimized))
        # Checking non-cds sequences
        else:
            self.optimized = self.original_sequence
            self.conditionals["cut_by"] = cutsite_list(part_type,self.optimized)
        # Apply to make new genbank file
        self.optimized_genbank = replace_genbank_sequence(self.original_genbank, self.optimized)
        # Check sequences
        if self.buildable():
            print(self.gene_id + " is clear of enzymes")
        else:
            print(self.find_enzyme() + " found in " + self.gene_id)
        if self.part_type == "cds" and not genbank_dictionary == {}:
            if transl_checker(self.optimized[:-3], self.genbank_dictionary["translation"]):
                print(self.gene_id + " translation is correct")
            else: 
                print(self.gene_id + " INCORRECT TRANSLATION")
                sys.exit()
        # Add cloning enzymes
        if len(self.optimized) < 300:
            self.cloning_enzyme = "BtgZI"
        else:
            self.cloning_enzyme = "BbsI"
        self.retrieval_enzyme = "BsaI"


    # Functions
    def buildable(self):
        enzyme = str(self.find_enzyme())
        if enzyme == "Clear":
            return True
        else:
            return False
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
            if single_finder(enzyme[1],seq) and single_finder(reverse_complement(enzyme[1]),seq):
                return enzyme[0]
            else:
                return "Clear"
    def dictionary_to_json_genbank(self,template, dictionary):
        for key, value in dictionary.items():
            template["genbank"][key] = dictionary[key]
        return template
    def json_write(self): 
        if self.buildable():
            path = "./../stage/{}".format(self.gene_id)
            os.makedirs(path)
            self.template_json["gene_id"] = self.gene_id
            self.template_json["author"]["name"] = self.author_name
            self.template_json["author"]["email"] = self.author_email
            self.template_json["author"]["affiliation"] = self.author_affiliation
            self.template_json["author"]["orcid"] = self.author_orcid
            self.template_json["info"]["documentation"]["gene_name"] = self.gene_name
            self.template_json["info"]["documentation"]["description"] = self.description
            self.template_json["info"]["documentation"]["project_description"] = self.project_description
            self.template_json["info"]["documentation"]["database_links"] = self.database_links
            self.template_json["info"]["gene_metadata"]["cloning"]["part_type"] = self.part_type
            self.template_json["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] = self.cloning_enzyme
            self.template_json["info"]["gene_metadata"]["cloning"]["retrieval_enzyme"] = self.retrieval_enzyme
            self.template_json["info"]["gene_metadata"]["cloning"]["target_organism"]["organism_name"] = self.target_organism
            self.template_json["info"]["gene_metadata"]["safety"] = self.safety
            self.template_json["info"]["gene_metadata"]["collection_id"] = self.collection_id
            self.template_json["info"]["gene_metadata"]["conditions"] = self.conditionals
            self.template_json["dates"]["submitted"] = self.submission_timestamp
            # Write taxid, target_organism
            self.template_json["sequence"]["original_sequence"] = self.original_sequence
            self.template_json["sequence"]["optimized_sequence"] = self.optimized
            self.template_json["tags"] = self.tags
            self.template_json["hashtags"] = self.hashtags
            if self.genbank_dictionary:
                self.dictionary_to_json_genbank(self.template_json, self.genbank_dictionary)
            with open("{}/{}.json".format(path,self.gene_id),"w+") as json_file:
                json.dump(self.template_json,json_file,indent=2)
            with open("{}/{}.gb".format(path,self.gene_id),"w+") as genbank_single:
                genbank_single.write(self.optimized_genbank)
        else: 
            print(self.gene_id + " NOT BUILDABLE")



### =============
### =============
### ORDER MANAGER
### =============
### =============
def freegenes_order():
    stage = config["STAGE_PATH"]
    
    small_seq_ids = []
    small_seqs = []
    large_seq_ids = []
    large_seqs = []
    
    def order_fragment_genes():
        for file in glob.glob(stage + "*/*.json"):
            with open(file,"r") as json_file:
                data = json.load(json_file)
                part_type = data["info"]["gene_metadata"]["cloning"]["part_type"]
                part_type = part_type.lower()
                part_type = part_type.replace(" ", "_")
                ortho_pair = ("","")
                if not part_type == "vector":
                    if not data["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] == "BtgZI":
                        ortho_pair = get_pair(data["gene_id"])
                fragments = FG_standard_fragment(data["sequence"]["optimized_sequence"],part_type,data["info"]["gene_metadata"]["cloning"]["cloning_enzyme"], ortho_pair)
                print(fragments)
                print(data["gene_id"])
                gene_id = data["gene_id"]
                for index,frag in enumerate(fragments):
                    fragment_name = gene_id + "_" + str(index + 1)
                    data["sequence"]["fragment_sequences"][fragment_name] = frag
                path = "{}{}".format(stage,gene_id)
                with open("{}/{}.json".format(path,gene_id),"w+") as json_file:
                    json.dump(data,json_file,indent=2)
    
    def order_write_link():
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
            if repeat_check(row["Sequence"] + small_row["Sequence"]):
                joined_seq = row["Sequence"] + small_row["Sequence"]
                joined_ids.append(row["Gene ID"])
                joined_seqs.append(joined_seq)
                fragment_names.append(row["Gene ID"] + "_link_" + small_row["Gene ID"])
                joined_ids.append(small_row["Gene ID"])
                joined_seqs.append(joined_seq)
                fragment_names.append(row["Gene ID"] + "_link_" + small_row["Gene ID"])
                small_counter += 1
            else:
                print(row["Gene ID"] + " AND " + small_row["Gene ID"])
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
    
    
    def order_twist_order():
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
    
    def order_replace_bad_sequence(gene_id):
        table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)
        json_data = Json_load(stage + "{}/{}.json".format(gene_id,gene_id))
        seq_to_replace = input("Sequence to replace? ")
        new_seq = input("New sequence? ")
        old_sequence = json_data["sequence"]["optimized_sequence"]
        new_sequence = fix_sequence(table,gene_id,old_sequence.replace(seq_to_replace, new_seq),0,[])
        if not Seq(old_sequence, IUPAC.unambiguous_dna).translate() == Seq(new_sequence, IUPAC.unambiguous_dna).translate():
            print("Bad translation, try again")
            replace_bad_sequence(gene_id)
        else:
            with open(stage + "{}/{}.gb".format(gene_id,gene_id),"r") as genbank_single:
                genbank_current = genbank_single.read()
            genbank_fixed = replace_genbank_sequence(genbank_current, new_sequence)
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
                order_fragment_to_order()
    
    def order_reoptimize_fragment(gene_id):
        json_data = Json_load(stage + "{}/{}.json".format(gene_id,gene_id))
        table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)
        translation = Seq(json_data["sequence"]["optimized_sequence"], IUPAC.unambiguous_dna).translate()[:-1]
        #translation = json_data["genbank"]["translation"]
        new_sequence = fix_sequence(table,gene_id,codon.optimize_protein(table, translation) + "TGA")
        if not Seq(json_data["sequence"]["optimized_sequence"], IUPAC.unambiguous_dna).translate() == Seq(new_sequence, IUPAC.unambiguous_dna).translate():
            print("Bad translation, try again")
            order_reoptimize_fragment(gene_id)
        else:
            with open(stage + "{}/{}.gb".format(gene_id,gene_id),"r") as genbank_single:
                genbank_current = genbank_single.read()
            genbank_fixed = replace_genbank_sequence(genbank_current, new_sequence)
            with open(stage + "{}/{}.gb".format(gene_id,gene_id),"w+") as genbank_single:
                genbank_single.write(genbank_fixed)
            json_data["sequence"]["optimized_sequence"] = new_sequence
            with open(stage + "{}/{}.json".format(gene_id,gene_id),"w+") as json_file:
                json.dump(json_data,json_file,indent=2)
            print("Wrote new sequence for " + gene_id)
    
    def order_reset_fragment_stage():
        for file in glob.glob(stage + "*/*.json"):
            with open(file,"r") as json_file:
                data = json.load(json_file)
                data["sequence"]["fragment_sequences"] = {}
                gene_id = data["gene_id"]
                path = "{}{}".format(stage,gene_id)
                with open("{}/{}.json".format(path,gene_id),"w+") as json_file:
                    json.dump(data,json_file,indent=2)
        print("Fragments on stage cleared.")
    
    
    def order_id_reset():
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
    


    
    def order_fragment_to_order():
        print("Recreating database")
        order_reset_fragment_stage()
        order_fragment_genes()
        order_write_link()
        order_twist_order()

    def order_manager():
        print("\n")
        print("=== FreeGenes Order manager ===")
        options = ("Frag -> Order", "Set ID starting point", "Fragment genes in stage", "Write linkers to stage", "Create Twist submission spreadsheet", "Change part of sequence", "Reoptimize single fragment", "Reset fragment stage", "Clear ortho pairs","Bad seq removal", "Reurn to FG manager", "Exit")
        choice = option_list(options)

        if choice == "Frag -> Order":
            order_fragment_to_order()
        elif choice == "Set ID starting point":
            order_id_reset()
        elif choice == "Fragment genes in stage":
            order_fragment_genes()
        elif choice == "Write linkers to stage":
            order_write_link()
        elif choice == "Create Twist submission spreadsheet":
            order_twist_order()
        elif choice == "Change part of sequence":
            replace_bad_sequence(input("gene_id : "))
        elif choice == "Reoptimize single fragment":
            reoptimize_fragment(input("gene_id : "))
        elif choice == "Reset fragment stage":
            order_reset_fragment_stage()
        elif choice == "Bad seq removal":
            order_bad_seq()
        elif choice == "Reurn to FG manager":
            return freegenes_manager()
        elif choice == "Exit":
            sys.exit()
        elif choice == "Clear ortho pairs":
            clear_pairs()

        print("Returning to Order manager")
        return order_manager()
    order_manager()


### ===============
### ===============
### GENOME DIGESTER
### ===============
### ===============

def freegenes_digest():
    ## ===============
    ## SETUP VARIABLES
    ## ===============
    date = datetime.date.today().strftime("%d") + "-" + datetime.date.today().strftime("%B")[:3].upper() + "-" + datetime.date.today().strftime("%Y")
    stage = FreeGenes_configuration()["STAGE_PATH"]


    NextCollection = NextCollection()
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

    important_tags = important_tags

    for file in glob.glob("./../pipeline/template.json"):
            with open(file,"r") as template_json:
                template = json.load(template_json)

    ## ================
    ## Digest functions
    ## ================
    # Process a string of important values into a dictionary
    def dictionary_builder(tag_list, string):
        return {key: value for (key, value) in map(lambda tag: [tag, ''.join(re.findall(r'/{}=\"([A-Za-z0-9:_./-_\s-]+)\"'.format(tag),string))], tag_list)}

    # Build dictionary into genbank compatible format
    def dictionary_to_genbank(dictionary):
        value_list = []
        for key, value in dictionary.items():
            if type(value) == type(""):
                value_list.append(str(key + "=" + '"' + value + '"'))
        return value_list

    def genbank_multiline(genbank_list):
        multiline= ''
        for item in genbank_list:
            split_list = textwrap.wrap(item, width=58)
            for item in split_list:
                multiline = multiline + "                     " + "/" + item + "\n"
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
    path = FreeGenes_configuration()["PROJECT_PATH"]

    print("Please choose a genome file")
    genome = file_list(path + "genome/genome_sequences/*.gb")
    print("Please choose configuration file")
    digest_configuration = file_list(path + "genome/digestion_configuration/*.yaml")
    config = load_configuration(digest_configuration)
    transl_table = '                     /transl_table='+str(config["transl_table"])


    ## =====================================
    ## Check for protein files or gene lists
    ## =====================================
    if config["protein_file"]:
        protein_fasta = file_list(path + "genome/protein_lists/*.fasta")
        protein_dictionary = fasta_refseq_dictionary(protein_fasta)

    if config["gene_list"]:
        essential_list = file_list(path + "genome/gene_lists/*.txt")
        essential_list = genelist_to_list(essential_list)

    ## ========================
    ## Process Genbank into CSV
    ## ========================
    df = genbank_to_csv(genome)


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
        multiline = (genbank_multiline(dictionary_to_genbank(data)) + transl_table + "\n" + "                     /codon_start=1" + "\n" + genbank_multiline(references))
        if not data["gene"] == "":
            definition = data["gene"]
        else:
            definition = data["locus_tag"]

        # Get gene_id and collection id
        gene_id = NextID()

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
        genbank_file += suffix_genbank(sequence)
        # Create object
        if write:
            freegene = FreeGene(gene_id, NextCollection, date, config["author"], config["email"], "BioBricks Foundation", "NA", definition, config["description"], links, "CDS", data["Source"], config["target_organism"], config["safety"], genbank_file, Json_load("./configuration/template.json"), config["optimization_table"], data, config["tags"])
            freegene.json_write()
        else:
            print('skipping ' + definition)

    return freegenes_manager()

## ========
## ========
## RETRIEVE
## ========
## ========
def freegenes_retrieve():

    ## ==========
    ## Parameters
    ## ==========
    
    fragment_condition = pd.DataFrame(columns=["Author", "Email", "Gene", "Error", "Genbank", "FG Pass", "Twist Pass", "Sent note", "Order confirmation"])

    def add_condition(author,email,gene,error,genbank,fg_pass,twist_pass=None, sent_note=None):
        new_fragment_condition = fragment_condition.append({"Author":author, "Email":email,"Gene":gene,"Error":error, "Genbank":genbank, "FG Pass":fg_pass, "Twist Pass":twist_pass, "Sent note":sent_note}, ignore_index=True)
        return new_fragment_condition

    counter = 0
    base_pairs = 0
    write = input("Write files to system? (ENTER to continue, 'WRITE' to write, 'UPDATE' for dashboard update) : ")

    template = Json_load("./configuration/template.json")
    config = FreeGenes_configuration()
    SINGLE_SUBMISSION_ID = config["SINGLE_SUBMISSION_ID"]
    BULK_SUBMISSION_ID = config["BULK_SUBMISSION_ID"]
    PREVIOUS_CSV_SINGLE = config["PREVIOUS_CSV_SINGLE"]
    PREVIOUS_CSV_BULK = config["PREVIOUS_CSV_BULK"]
    UPDATE_CSV_BULK = config["UPDATE_CSV_BULK"]
    UPDATE_CSV_SINGLE = config["UPDATE_CSV_SINGLE"]
    ## ======
    ## Script
    ## ======

    # Bulk data
    current_data_bulk = fill_strip(extract_google_form(BULK_SUBMISSION_ID))
    previous_bulk = fill_strip(pd.read_csv(PREVIOUS_CSV_BULK))
    bulk_data = uniq_data(current_data_bulk, previous_bulk)
    for index, row in bulk_data.iterrows():
        collection_id = NextCollection()
        csv_url = row['CSV file']
        csv_data = csvtext_to_pandas(get_wufoo_textfile(csv_url)).fillna('')
        project_description = row['project_description']
        for index_csv, row_csv in csv_data.iterrows():
            # Setup import into object
            gene_id = NextID(counter) # New ID
            zip_file = zipfile.ZipFile(io.BytesIO(requests.get(row['Genbank files']).content)) # Download zip file
            name_list = [ x for x in zip_file.namelist() if "__MACOSX" not in x ]
            genbank_file = ""
            for name in name_list:
                if row_csv['Exact genbank name in zip file'] == name:
                    with zip_file.open(name) as myfile:
                        genbank_file="\n".join(str(myfile.read(), 'utf-8').splitlines()) # Format genbank file all nice
            ## IMPLEMENT CONDITIONAlS
            conditionals = {}
            try:
                tags = row_csv['tags']
            except:
                tags = []
            try:
                hashtags = row_csv['Hashtags'].split("#")
            except:
                hashtags = []
            try:
                optimize = row_csv["optimize"]
            except:
                optimize = False
            try:
                # Import into object
                freegene = FreeGene(gene_id, collection_id, row['date'], row['Name'], row['Email'], row['Affiliation'], row['ORCID'], row_csv['Gene name'], row_csv['Description'], row_csv['Links'], row_csv['Part type'], row_csv['Source organism'], row_csv['Target organism'], row_csv['Safety information'], genbank_file, template, optimize, project_description=project_description, tags=tags, hashtags=hashtags, conditionals=conditionals)
                if write == "WRITE":
                    freegene.json_write()
                else:
                    counter += 1
                base_pairs = base_pairs + len(freegene.optimized)
                error = "No data errors."
                FG_pass = retrieve_checker(freegene.optimized)
            except Exception as e:
                error = e
                FG_pass = "Failed entry."
            fragment_condition = add_condition(row['Name'],row['Email'],row_csv['Gene name'],error,genbank_file,FG_pass)
    if write == "WRITE":
        current_data_bulk.to_csv(PREVIOUS_CSV_BULK, index=False)
    if write == "UPDATE":
        bulk_data.to_csv(UPDATE_CSV_BULK, index=False)

    # Single submissions
    collection_id = NextCollection()
    current_data_single = fill_strip(extract_google_form(SINGLE_SUBMISSION_ID))
    previous_single = fill_strip(pd.read_csv(PREVIOUS_CSV_SINGLE))
    single_data = uniq_data(current_data_single, previous_single)
    for index, row in single_data.iterrows():
        gene_id = NextID(counter) # New ID
        genbank_file = get_wufoo_textfile(row['genbank_file'])
        if row['optimize'] == "Yes":
            optimize = "custom_1"
        else:
            optimize = False
        try:
            freegene = FreeGene(gene_id, collection_id, row['Timestamp'], row['name'], row['email'], row['affiliation'], row['orcid'], row['gene_name'], row['description'], row['links'], row['part_type'], row['source_organism'],  row['target_organism'], row['safety_info'], genbank_file, template, optimize, project_description=row["project_description"], tags=row["tags"].split(","), hashtags=row["hashtags"].split("#"), conditionals={})
            if write == "WRITE":
                freegene.json_write()
            else: 
                counter += 1
            base_pairs = base_pairs + len(freegene.optimized)
            error = "No data errors."
            FG_pass = retrieve_checker(freegene.optimized)
        except Exception as e:
            error = e
            FG_pass = "Failed entry."
        fragment_condition = add_condition(row['name'],row['email'],row['gene_name'],error,genbank_file,FG_pass)
    if write == "WRITE":
        current_data_single.to_csv(PREVIOUS_CSV_SINGLE, index=False)
    if write == "UPDATE":
        single_data.to_csv(UPDATE_CSV_SINGLE, index=False)

    # Text output
    print("Total of " + str(counter) + " genes to synthesize")
    print("Numbering approximately " + str(base_pairs) + " base pairs")
    #gmail_password = getpass.getpass("Enter Email Password: ")
    #freegenes_sendmail(gmail_password, ['koeng101@gmail.com'], "Failed sequences", failed)


    # Only update dataframe with new entries
    fragment_condition.to_csv("./configuration/gene_status/fragment_condition.csv")

    # Update the dashboard ADD STUFF TO CONFIG FOR BIONET SYNTHESIS
    if write == "UPDATE":
        with open(config['BIONET_DASHBOARD'], 'w') as f:
            f.write(update_dashboard())
        print("Remember to push to git!")

    return freegenes_manager()

## ===========
## Stage email
## ===========

def email_manager():
    print("\n")
    print("=== FreeGenes Email Manager ===")
    options = ("Send update from latest retrieve","Send Twist synthesis update","Exit")
    choice = option_list(options)



def email_update():
    status = pd.read_csv("./configuration/gene_status/fragment_condition.csv")
    print(status)


def write_email(author,email,genes,errors,fg_pass,twist_pass):
    hello = '''

    '''
    return hello

##
## Dashboard
##

def update_dashboard():
    def process_strings():
        df_single = pd.read_csv(config["UPDATE_CSV_SINGLE"])[['Timestamp', 'name', 'gene_name', 'description', 'links']].copy()
        single_submission_str = ''
        for submission in df_single.values.tolist():
            for item in submission:
                single_submission_str += item.replace('\n', ' ') + ' | '
            single_submission_str += "\n"

        df_bulk = pd.read_csv(config["UPDATE_CSV_BULK"])[['date', 'Name', 'project_description']].copy()
        bulk_submission_str = ''
        for submission in df_bulk.values.tolist():
            for item in submission:
                bulk_submission_str += item.replace('\n', ' ') + ' | '
            bulk_submission_str += "\n"

        return [single_submission_str, bulk_submission_str]

    submission_list = process_strings()
    return '''
# Simple dashboard

All single requests can be found [here](https://docs.google.com/spreadsheets/d/1j5Gc7KEfRlPCIaXMGjDhgQDfSOVx7tnbss9AksrHhzk/edit?usp=sharing). 

All bulk requests can be found [here](https://docs.google.com/spreadsheets/d/1qgNt3h63--o7qlhTdkUqGpLqVOGizUaY5dMG7VdCwHA/edit?usp=sharing).

## New single requests:

Date | Author | Gene Name | Description | Links
| --- | --- | --- | --- | --- |
{}
## New bulk requests: 

Date | Author | Project Description 
| --- | --- | --- |
{}
'''.format(submission_list[0], submission_list[1])

    


## =================
## =================
## FreeGenes manager
## =================
## =================

def freegenes_manager():
    print("\n")
    print("=== FreeGenes Manager ===")
    options = ("Retrieve from order sheet","Digest a genome","Order DNA","Email update","Exit")
    choice = option_list(options)

    if choice == "Retrieve from order sheet":
        return freegenes_retrieve()
    elif choice == "Digest a genome":
        return freegenes_digest()
    elif choice == "Order DNA":
        return freegenes_order()
    elif choice =="Email update":
        return email_update()
    elif choice == "Exit":
        sys.exit()
    else:
        return freegenes_manager()

