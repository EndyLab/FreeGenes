import collections
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


## ==========
## Parameters
## ==========

cut_sites= [
        ("BbsI", "GAAGAC"),
        ("BtgZI", "GCGATG"),
        ("BsaI", "GGTCTC"),
        ("BsmBI", "CGTCTC"),
        ("AarI", "CACCTGC"),
        ("BfuAI", "ACCTGC")]

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


## ==============
## Gene functions
## ==============

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATGC","TACG"))[::-1]

def translate(seq):
    return str(Seq(seq, generic_dna).translate(table=11))

def random_dna_sequence(length):
    return ''.join(np.random.choice(('A', 'C', 'T', 'G')) for _ in range(length))

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

def transl_checker(seq, aa):
    prot = translate(seq)
    if prot[-1] == "*":
        prot = prot[:-1]
    return prot == aa


## ===============
## Recode_sequence
## ===============

def recode_sequence(table, seq, rep, halfway=False):
    if seq.find(rep) < 0:
        return seq

    # This block contains code that picks a new codon.
    def list_del(del_list, position):
        del del_list[position]
        return del_list
    def fraction_list(number_list):
        return list(map(lambda x: x / sum(number_list), number_list))
    def new_codon_chooser(old_codon):
        # Add exception handling for when there are 0 cases of said codon. 
        recode_aa = str(Seq(old_codon, IUPAC.unambiguous_dna).translate())
        old_codon_location = list(table.loc[recode_aa].index).index(old_codon.upper())
        return ''.join(choice(list_del(list(table.loc[recode_aa].index),old_codon_location), 1, p=fraction_list(list_del(list(table.loc[recode_aa]['Number']),old_codon_location))))

    # This block contains code that produces the desired position information
    def choose_position(list_of_positions, halfway):
        if halfway:
            position = list_of_positions[int(len(list_of_positions)/2)]
        else:
            position = list_of_positions[0]
        return position
    def position_list(seq, rep):
        position = seq.find(rep)
        position -= position % 3
        list_of_positions = list(range(position, position + (len(rep) // 3 + 1) * 3, 3))
        return list_of_positions

    position_choice = choose_position(position_list(seq, rep), halfway)
    old_codon = seq[position_choice:position_choice+3]
    new_codon = new_codon_chooser(old_codon)
    print("{} -> {} at position {}-{}.".format(old_codon,new_codon,position_choice+1,position_choice+3))
    return seq[:position_choice] + new_codon + seq[position_choice+3:]


def FG_MoClo_ends(seq):
    return seq[:-6].join(list(map(lambda x: end_codons[str(Seq(x, IUPAC.unambiguous_dna).translate())], [seq[-6:-3],seq[-3:]])))



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
    return len(repeat_finder(seq)) < 20

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

## =====
## Fixer
## ===== 

def fix_sequence(table,gene_id,seq,fix_attempts=0):
    fix_attempts = fix_attempts + 1
    max_attempts = 1000
    if fix_attempts > max_attempts:
        print("{} could not be fixed. Please manually check this sequence.".format(gene_id))
        print("Terminating program")
        sys.exit()
    if not gc_in_range(seq):
        print("Fixing GC content in {}".format(gene_id))
        return fix_sequence(table, gene_id, codon.optimize_protein(table, translate(seq)))
    elif not homopolymer_checker(seq):
        print("Fixing homopolymer stretches in {}".format(gene_id))
        homopolymer_max, homopolymer_site = max_homopolymer(seq)
        homopolymer = seq[homopolymer_site:homopolymer_site+homopolymer_max]
        return fix_sequence(table, gene_id, recode_sequence(table, seq, homopolymer, halfway=True), fix_attempts)
    elif not repeat_check(seq):
        print("Fixing repeat regions in {}".format(gene_id))
        return fix_sequence(table, gene_id, recode_sequence(table, seq, repeat_finder(seq), halfway=True), fix_attempts)
    elif not cutsite_check(seq):
        print("Fixing enzyme cutsites in {}".format(gene_id))
        return fix_sequence(table, gene_id, recode_sequence(table, seq, dict((y, x) for x, y in cut_sites)[find_enzyme(seq)], halfway=False), fix_attempts)
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
            return fix_sequence(table,gene_id,seq)
    print("{} successfully checked.".format(gene_id))
    return seq

    
table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)
seq = "ATGGAGAAGAATCGTTCCGCTTTTCAGCAAAACCAGCAAGCCAGCAACCAGCCTTTCAATCAAGACCAGAACCAATATTATCAGGACCCGAACCAACAACAGTTTAATCAATCAGGTTTCGATCCGAATCAGCAGCAATTCAACCAGCCGGGTTTCGATCCTAACCAGCAATACTACCAGGATCCGAATCAGCAGCAATTCAATCAGGCCGGCTTCGACCAAAATCAACAGTACTACCAGGACCCGAATCAGCAGCAATTTAACCAACCGGGCTTCGACCCGAATCAACAGTACTACCAGGACCCAAACCAACAGCAGTTCAACCAGGCAGGCTTTGACCAAAATCAATACTACCAAGACCCGAACCAACAACAGTTCAATCAAAGCGGTTTCGATCAGAACCAATACTACCAAGATCCAAACCAGCAACAGTTCAACCAACCGTCCTTCGACCTGAATAATCAGCAGTTCAACCAGCCGGGTTTCAATCAAAGCCCTGCTTTCGAAATCACCCCGCAGGAACAAAAAGCTGAGCAGGAAATGTTCGGTGAGGAGCCGCCGCAAGTTGTTCGCGAAATTCACGAACTGCCATTTGAGAAAATTCGTAGCTTCCTGCAGAGCGACTTTGATAGCTACAATTTTCGTCTGAATTCCCTGAAATCTAAGTTGGACAATGCGCTGTATAGCCTGGATAAAACCATTCAAAACACCAACGAGAACACGGCGAACCTGGAGGCGATCCGTCACAACCTGGAACAGAAAATCCAAAATCAGAGCAAACAGCTGCGCACCAACTTTGACACGCAGAAACTGGATGACAAAATTAATGAGCTGGAAATTCGTATGCAGAAACTGACCCGTAACTTCGAAAGCCTGAGCGAATTGTCCAAGCATAACAGCTACCCGAACTACTACGAGAAACTGTTGCCGAACGGTGGCGATTCTATGACCAATGTGTTTGAGAAGGCACTGATGATGAATCTGCTGCGTACTACGCTGCCGCCACAGCCACAGGTTCAATACTACCCGCAGCCGTACCCGTATATTCGTCCGTACTATGACGAACCGATCTACGCCGGTTTCCGTCGTCGTGGCTACCGTGACGACTTTTATGAGTGA"


new_seq = check_sequence(table,"my gene",seq)
print(new_seq)

