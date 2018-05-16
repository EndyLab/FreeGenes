import numpy as np
import pandas as pd

pairs_file = ""
primer_pairs = pd.DataFrame()

def load_pairs(file = "./orthogonal_primers_pairs.csv"):
    global primer_pairs
    global pairs_file
    
    pairs_file = file
    primer_pairs = pd.read_csv(pairs_file, index_col=0)

def save_pairs():
    primer_pairs.to_csv(pairs_file)
    
def clear_pairs():
    primer_pairs['GeneID'] = np.nan
    save_pairs()
    
def get_pair(gene_id):
    index = primer_pairs['GeneID'].isnull().idxmax()
    primer_pairs.loc[index, 'GeneID'] = gene_id
    save_pairs()
    return (
        primer_pairs.loc[index, 'Forward'],
        primer_pairs.loc[index, 'Reverse']
    )