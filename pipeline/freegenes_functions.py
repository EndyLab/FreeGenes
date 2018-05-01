import collections
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
import codon
from Bio.Alphabet import IUPAC
import yaml
from numpy.random import choice


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

def NextID():
    counter = len(glob.glob("./../stage/" + "*")) 
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

def fix_sequence(table,gene_id,seq,fix_attempts=0,change_list=[]):
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
    elif not gc_50bp_acceptable_range(seq):
        print("Fixing 50bp of outlier GC in {}".format(gene_id))
        gc_outlier = gc_range_sequence(seq)
        new_seq = recode_region(table,gene_id,seq,gc_outlier)
        return fix_sequence(table, gene_id, new_seq , fix_attempts, change_list)
    elif not gc_in_range(seq):
        print("Fixing GC content in {}".format(gene_id))
        return fix_sequence(table, gene_id, codon.optimize_protein(table, translate(seq)))
    elif not univerisal_primer_checker(seq):
        print("Fixing internal FG primer regions in {}".format(gene_id))
        return fix_sequence(table, gene_id, universal_primer_finder(seq), fix_attempts, change_list)
    else:
        return check_sequence(table,gene_id,seq)


## =======
## Checker
## =======

def check_sequence(table,gene_id,seq):
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
    sequence_checks = [is_seq, is_triplet, has_start, has_no_internal_stops, MoClo_ends, homopolymer_checker, gc_in_range, cutsite_check, repeat_check]
    for test in sequence_checks:
        if not test(seq):
            print("{} failed on {}".format(gene_id, test.__name__))
            input("PRESS ENTER TO FIX THIS GENE!")
            return fix_sequence(table,gene_id,seq)
    print("{} CHECKED.".format(gene_id))
    return seq

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

## ================================
## Twist fix gene by reoptimization
## ================================
def replace_genbank_sequence(original_genbank, optimized_sequence):
    return original_genbank.replace(''.join(re.findall(r'(ORIGIN[A-Za-z0-9:_./-_\s-]+//)', original_genbank)), suffix_genbank(optimized_sequence))

def reoptimize(gene_id):
    json_data = Json_load("./../data/{}/{}.json".format(gene_id,gene_id))
    new_sequence = final_fixer(json_data["genbank"]["translation"])
    with open("./../data/{}/{}.gb".format(gene_id,gene_id),"r") as genbank_single:
        genbank_current = genbank_single.read()
    genbank_fixed = replace_genbank_sequence(genbank_current, new_sequence)
    with open("./../data/{}/{}.gb".format(gene_id,gene_id),"w+") as genbank_single:
        genbank_single.write(genbank_fixed)
    json_data["sequence"]["optimized"] = new_sequence
    with open("./../data/{}/{}.json".format(gene_id,gene_id),"w+") as json_file:
        json.dump(json_data,json_file,indent=2)
    print("Wrote new sequence for " + gene_id)
    print("Remember to refragment, order, and make Twist submission")


def final_fixer(translation):
    optimize = "custom_1"
    #fix_sequence = codon.optimize_protein(codon.load_codon_table(taxonomy_id=optimize, custom=True), translation) + "TAA"
    optimized = optimize_gene("gene", fix_sequence, "83333")#self.target_organism)
    return optimized



## =======
## Classes
## =======
class FreeGene:
    """FreeGene class"""
    def __init__(self, gene_id, collection_id, timestamp, author_name, author_email, author_affiliation, author_orcid, gene_name, description, database_links, part_type, source_organism, target_organism, safety, genbank_file, template_json, optimize, genbank_dictionary={}, tags=[]):
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
        self.genbank_dictionary = genbank_dictionary
        self.tags = tags
        self.original_genbank = genbank_file
        self.original_sequence = str(SeqIO.read(StringIO(genbank_file), "genbank").seq)
        # Sequence modifications
        print("Running {}".format(self.gene_name))
        if self.part_type == "cds":
            ####### FIX: ONLY DOES CODON CHANGING ACCORDING TO THE E COLI TABLE!!
            if optimize == False:
                optimized_sequence = self.original_sequence
                optimization_table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)
            else:
                print("\n\nOptimizing {}. ( {} )".format(self.gene_id, self.gene_name))
                optimization_table = codon.load_codon_table(taxonomy_id=optimize, custom=True)
                optimized_sequence = codon.optimize_protein(optimization_table, genbank_dictionary["translation"]) + "TAA"
            print(optimized_sequence)
            self.optimized = fix_sequence(optimization_table,self.gene_id,optimized_sequence,0,[])
        else:
            self.optimized = self.original_sequence
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
            self.template_json["info"]["documentation"]["database_links"] = self.database_links
            self.template_json["info"]["gene_metadata"]["cloning"]["part_type"] = self.part_type
            self.template_json["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] = self.cloning_enzyme
            self.template_json["info"]["gene_metadata"]["cloning"]["retrieval_enzyme"] = self.retrieval_enzyme
            self.template_json["info"]["gene_metadata"]["cloning"]["target_organism"]["organism_name"] = self.target_organism
            self.template_json["info"]["gene_metadata"]["safety"] = self.safety
            self.template_json["info"]["gene_metadata"]["collection_id"] = self.collection_id
            self.template_json["dates"]["submitted"] = self.submission_timestamp
            # Write taxid, target_organism
            self.template_json["sequence"]["original_sequence"] = self.original_sequence
            self.template_json["sequence"]["optimized_sequence"] = self.optimized
            self.template_json["tags"] = self.tags
            if self.genbank_dictionary:
                self.dictionary_to_json_genbank(self.template_json, self.genbank_dictionary)
            with open("{}/{}.json".format(path,self.gene_id),"w+") as json_file:
                json.dump(self.template_json,json_file,indent=2)
            with open("{}/{}.gb".format(path,self.gene_id),"w+") as genbank_single:
                genbank_single.write(self.optimized_genbank)
        else: 
            print(self.gene_id + " NOT BUILDABLE")


