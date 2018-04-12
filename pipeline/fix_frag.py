import codon
import freegenes_functions as ff
def FG_MoClo_gene(sequence,table):
    def pick_codon(aa):
        return codon.optimize_protein(table, aa) 

    def recode_sequence(seq, rep):
        pos = seq.find(rep)
        if pos < 0:
            return seq
        pos -= pos % 3
        for i in range(pos, pos + (len(rep) // 3 + 1) * 3, 3):
            codon = seq[i:i+3]
            choices = ec_codon_usage_10plus.ix[ec_codons.ix[codon]]
            choices = choices[choices.index != codon]
            if choices.shape[0] > 0:
                newcodon = choices.iloc[(choices.Fraction.cumsum() / choices.Fraction.cumsum().max() < np.random.rand()).sum()].name # Stochastically allocate codon
    #             newcodon = choices[choices.index != codon].Fraction.idxmax() # Deterministically allocate codons by using the most frequence one
                break

        eprint("{} -> {}".format(codon, newcodon))

        return seq[:i] + newcodon + seq[i+3:]

    def remove_cutsites(name, seq):
        changes = 0
        for enzyme, cut in cut_sites + [(e, reverse_complement(c)) for e, c in cut_sites]:
            while cut in seq:
                eprint("{} cuts {} ({})".format(enzyme, name, cut))
                changes += 1
                seq = recode_sequence(seq, cut)
        return seq, changes

    for i, gene in design_genes.iterrows():
        eprint("Fixing gene {}".format(gene['Gene']))

        sequence = gene['Sequence']
        for j in range(1000): # try 1000 times to get a working sequence
            sequence, changed = remove_cutsites(gene['Gene'], sequence)

            homopolymer_length, homopolymer_pos = max_homopolymer(sequence)
            while homopolymer_length >= 6:
                changed = True
                sequence = recode_sequence(sequence, sequence[homopolymer_pos:homopolymer_pos+homopolymer_length])
                homopolymer_length, homopolymer_pos = max_homopolymer(sequence)

            if not gc_in_range(sequence):
                eprint ("GC out of range")

            if not changed:
                break
        else:
            eprint("Warning: {} could not be optimized".format(gene['Gene']))
        design_genes.ix[i, 'Sequence'] = sequence

    # Force final codons to our standard set
    force_codons = {
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
        'G': 'GGC'
    }

    sequence = sequence.str[:-6] + sequence.str[-6:-3].apply(lambda x: force_codons[ec_codons.ix[x]]) + "TGA"

    # Make sure we didn't introduce an RE site
    # Make sure we didn't introduce any restriction sites
    error = False
    for enzyme, cut in cut_sites + [(e, reverse_complement(c)) for e, c in cut_sites]:
        if design_genes['Sequence'].str.contains(cut).any():
            error = True
            eprint("Cannot normalize final codon without introducing restriction site {} ({})".format(enzyme, cut))

            for i, (g, s) in design_genes.ix[design_genes['Sequence'].str.contains(cut), ['Gene','Sequence']].iterrows():
                eprint(g, s[s.find(cut):])
            eprint()

    if error:
        again = input("Try again? ( y or n )").upper()


##
##
##


def fragment_gene(sequence,entry_type):
    ## ==========================================================
    ## Configurations
    ## ==========================================================

    random_sequence = { # For Twist multigene sets
            "seq": "ATACACCGAC"
            }
    synthesis_configuration = {
        "max_length" : 1500,
        "min_length" : 300
    }
    pcr_configuration = {
        "max_length" : 5000
    }
    standard_flanks = {
        "prefix" : "GGAG",
        "suffix" : "CGCT"
    }
    enzyme_configuration = {
        "AarI" : {
            "prefix" : "CACCTGCCCTA",
            "suffix" : "ACTCGCAGGTG"
            },
        "BbsI" : {
             "prefix" : "GAAGACTA",
             "suffix" : "ACGTCTTC"
             },
        "BfuAI" : {
            "prefix" : "ACCTGCCCTA",
            "suffix" : "ACTCGCAGGT"
            },
        "BsaI" : {
            "prefix" : "GGTCTCA",
            "suffix" : "CGAGACC"
            },
        "BsmBI" : {
            "prefix" : "CGTCTCA",
            "suffix" : "CGAGACG"
            },
        "BtgZI" : {
            "prefix" : "GCGATGCATGCTTGCA",
            "suffix" : "CCTCTGAATACATCGC"
            },
        "SapI" : {
            "prefix" : "GCTCTTCA",
            "suffix" : "GCTCTTCA"
            },
        "None" : {
            "prefix" : "",
            "suffix" : ""
            },
        "BsaI-methylated_prefix" : {
            "prefix" : "CCGGTCTCA",
            "suffix" : "CGAGACC"
            },
        "BsaI-methylated_suffix" : {
            "prefix" : "GGTCTCA",
            "suffix" : "CGAGACCGG"
            },

    }
    # Establish part types
    part_type = dict(
        cds = {
            # FreeGenes MoClo CDS definition
            "cloning_enzyme": "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme": "BsaI-methylated_suffix",
            "prefix": "A",
            "suffix": "AGAGCTT"
        },
        eukaryotic_promoter = {
            # FreeGenes MoClo Eukaryotic Promoter definition
            "cloning_enzyme" : "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme" : "BsaI",
            "prefix" : "GGAG",
            "suffix" : "AATG"
        },
        prokaryotic_promoter = {
            # FreeGenes MoClo Prokaryotic Promoter definition
            "cloning_enzyme" : "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme" : "BsaI",
            "prefix" : "GGAG",
            "suffix" : "TACT"
        },
        rbs = {
            # FreeGenes MoClo RBS definition
            "cloning_enzyme" : "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme" : "BsaI",
            "prefix" : "TACT",
            "suffix" : "AATG"
        },
        terminator = {
            # FreeGenes MoClo Prokaryotic Promoter definition
            "cloning_enzyme" : "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme" : "BsaI",
            "prefix" : "GCTT",
            "suffix" : "CGCT"
        },
        operon = {
            # FreeGenes MoClo Operon definition
            "cloning_enzyme" : "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme" : "BsaI",
            "prefix" : "GGAG",
            "suffix" : "CGCT"
        }
    )
    ## ================================
    ## Defined find_enzyme for vectors
    ## ================================
    def find_enzyme(seq):
        seq = str(seq)
        def reverse_complement(seq):
            return seq.translate(str.maketrans("ATGC","TACG"))[::-1]

        cut_sites = [
            ("BbsI", "GAAGAC"),
            ("BtgZI", "GCGATG"),
            ("BsaI", "GGTCTC"),
            ("BsmBI", "CGTCTC"),
            ("AarI", "CACCTGC"),
            ("BfuAI", "ACCTGC")]
        def single_finder(enzyme_cut,sequence):
            if enzyme_cut in sequence:
                return True
            else:
                return False
        for enzyme in cut_sites:
            if not single_finder(enzyme[1],seq) and not single_finder(reverse_complement(enzyme[1]),seq):
                return enzyme[0]

    ## ==========================================================
    ## Add Retrieval prefix and suffix
    ## ==========================================================
    if entry_type not in part_type.keys() and "vector" not in entry_type:
        print("not a valid type")
        return
    if "vector" not in entry_type:
        part = part_type[entry_type]
        retrieval_enzyme = part["retrieval_enzyme"]
        full_sequence = enzyme_configuration[retrieval_enzyme]["prefix"] + random_dna_sequence(1) + part["prefix"] + sequence + part["suffix"] + random_dna_sequence(1) + enzyme_configuration[retrieval_enzyme]["suffix"] # Have a programmatic way to condense prefix / suffix rather than just defining them
    else:
        full_sequence = sequence.upper()
    ## =======================================
    ## Add Cloning prefix and suffix
    ## =======================================
    if "vector" in entry_type:
        seq = full_sequence + full_sequence[:4]
        cloning_enzyme = find_enzyme(seq)
    else:
        seq =  full_sequence
        cloning_enzyme = part["cloning_enzyme"]
    num_frags = len(seq) // synthesis_configuration["max_length"] + 1
    frag_len = len(seq) // num_frags

    frags = []
    if len(seq) < 300:
        small_cloning_enzyme = part["small_cloning_enzyme"]
        frag = enzyme_configuration[small_cloning_enzyme]["prefix"] + standard_flanks["prefix"] + seq + random_sequence["seq"] + standard_flanks["suffix"] + enzyme_configuration[small_cloning_enzyme]["suffix"]
        frags.append(frag)
    else:
        for i in range(num_frags):
            frag = seq[max(0, i * frag_len - 2):min((i+1) * frag_len + 2,len(seq))]
            frag = enzyme_configuration[cloning_enzyme]["prefix"] + standard_flanks["prefix"] + frag + standard_flanks["suffix"] + enzyme_configuration[cloning_enzyme]["suffix"]
            frags.append(frag)
    return frags

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

def find_enzyme(seq):
    seq = str(seq)
    cut_sites = [
        ("BbsI", "GAAGAC"),
        ("BtgZI", "GCGATG"),
        ("BsaI", "GGTCTC"),
        ("BsmBI", "CGTCTC"),
        ("AarI", "CACCTGC"),
        ("BfuAI", "ACCTGC")]
    for enzyme in cut_sites:
        if not sequence_search(enzyme[1],seq):
            return enzyme[0]

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
            "prefix" : "GGAGGTCTCNA",
            "suffix" : "AGAGCTTNGAGACCGCT"
            },
        "eukaryotic_promoter" : {
            "prefix" : "GGAGGTCTCNGGAG",
            "suffix" : "AATGNGAGACCGCT"
            },
        "prokaryotic_promoter" : {
            "prefix" : "GGAGGTCTCNGGAG",
            "suffix" : "TACTNGAGACCGCT"
            },
        "rbs" : { 
            "prefix" : "GGAGGTCTCNTACT",
            "suffix" : "AATGNGAGACCGCT"
            },
        "terminator" : {
            "prefix" : "GGAGGTCTCNGCTT",
            "suffix" : "CGCTNGAGACCGCT"
            },
        "operon" : { 
            "prefix" : "GGAGGTCTCNGGAG",
            "suffix" : "CGCTNGAGACCGCT"
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

def FG_standard_fragment(seq, part_type, cloning_enzyme):
    random_forward = "CATGCTTGCA"
    random_reverse = "GCTCTGAATA"
    cloning_sites = enzymes[cloning_enzyme]["seq"] + ("N" * enzymes[cloning_enzyme]["jump"])
    cloning_prefix = cloning_sites.replace("N" * cloning_sites.count("N"), random_forward[:cloning_sites.count("N")]) 
    cloning_suffix = reverse_complement(cloning_sites).replace("N" * cloning_sites.count("N"), random_reverse[:cloning_sites.count("N")])
    return fragmenter(part_type_preparer(part_type, seq), cloning_prefix, cloning_suffix)


    

