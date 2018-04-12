#!/usr/bin/python

'''

 

 

Preamble. Codon use frequency for three different coding schemes (DNA 2.0 schemes for E. coli). Scheme A=1st number, B=second, C= third. I usually use scheme B.

 

 

'''

# Stuff. python endy_dna20_codons.py

import random

# Set vfasta to protein sequences to be codon optimized.

vfasta="C:\\Users\\cmerryma\\Desktop\\Mgen proteins.FASTA.txt"

# Set vtable_name to A,B, or C to select frequency table.

vtable_name="B"

# Set filename of codon optimized nt output.

vout="C:\\Users\\cmerryma\\Desktop\\Mgen_optimized.fas"

 

 

# Compact but moderately readable form of frequency tables.

# For corresponding paper, GooScho "Design Parameters to Control Synthetic Gene Expression in Escherichia coli"

v1=""

v1=v1+'A             GCA~0.15~0.24~0~~GCC~0.18~0.19~0~~GCG~0.51~0.44~0.73~~GCT~0.15~0.12~0.26||'

v1=v1+'C              TGC~0.65~0.58~1~~TGT~0.35~0.42~0||'

v1=v1+'D             GAC~0.59~0.54~0.61~~GAT~0.4~0.46~0.38||'

v1=v1+'E              GAA~0.41~0.43~0.36~~GAG~0.59~0.57~0.63||'

v1=v1+'F              TTC~0.58~0.55~0.5~~TTT~0.42~0.45~0.5||'

v1=v1+'G             GGA~0.01~0~0~~GGC~0.41~0.39~0.43~~GGG~0.01~0~0~~GGT~0.57~0.6~0.56||'

v1=v1+'H             CAC~0.58~0.62~0.5~~CAT~0.42~0.38~0.5||'

v1=v1+'I               ATA~0~0~0~~ATC~0.46~0.51~0.77~~ATT~0.53~0.49~0.22||'

v1=v1+'K              AAA~0.44~0.51~0.35~~AAG~0.56~0.49~0.64||'

v1=v1+'L CTA~0~0~0~~CTC~0.05~0.03~0.18~~CTG~0.61~0.78~0.27~~CTT~0.04~0.02~0~~TTA~0.06~0.03~0.11~~TTG~0.24~0.14~0.43||'

v1=v1+'M            ATG~1~1~1||'

v1=v1+'N             AAC~0.56~0.53~0.64~~AAT~0.44~0.47~0.35||'

v1=v1+'P              CCA~0.15~0.1~0.13~~CCC~0~0~0~~CCG~0.72~0.81~0.86~~CCT~0.13~0.09~0||'

v1=v1+'Q             CAA~0.51~0.45~0.53~~CAG~0.48~0.55~0.46||'

v1=v1+'R AGA~0.09~0.03~0~~AGG~0.01~0~0~~CGA~0~0~0~~CGC~0.4~0.35~0.2~~CGG~0~0~0~~CGT~0.49~0.62~0.8||'

v1=v1+'S AGC~0.38~0.68~0.16~~AGT~0.01~0~0~~TCA~0.02~0.02~0~~TCC~0.25~0.13~0.58~~TCG~0.21~0.05~0.04~~TCT~0.12~0.1~0.2||'

v1=v1+'T              ACA~0~0~0~~ACC~0.63~0.57~0.91~~ACG~0.25~0.33~0~~ACT~0.11~0.1~0.08||'

v1=v1+'V             GTA~0.08~0.03~0~~GTC~0.22~0.28~0~~GTG~0.38~0.34~0.4~~GTT~0.31~0.35~0.59||'

v1=v1+'W            TGG~1~1~1||'

v1=v1+'Y              TAC~0.57~0.58~0.57~~TAT~0.42~0.42~0.42'

 

vtblA={}

vtblB={}

vtblC={}

 

# Amino acids separated by double pipe.

v2=v1.split("||")

for a,b in enumerate(v2):

                v3=b[0:b.find("\t")]

                v4=b[b.find("\t")+1:]

 

                # Codons separated by double squiggly.

                v5=v4.split("~~")

                for c,d in enumerate(v5):

                                v6=d[0:d.find("~")]

                                v7=d[d.find("~")+1:]

 

                                # Three Frequency sets separated by squiggly.

                                # Expand to 100 total choices of synonymous codons for each aa. To be randomized later.

                                # Eventually, we will pick one codon from the 100 as a protein is transformed amino acid by amino acid.

                                v8=v7.split("~")

                                for e,f in enumerate(v8):

                                                vA=[]

                                                v9=int(float(f)*100)

 

                                                # Skip codons that aren't used at all.

                                                if(v9==0):

                                                                continue

 

                                                # Or, generate its representation.

                                                else:

                                                                for g in range(0,v9):

                                                                                vA.append(v6)

 

                                                                # I'm hashing all 3 tables as I can foresee situations where that could be useful with additional code.

                                                                # Amino acid one letter code -> synonymous codon representation.

                                                                if(e==0):

                                                                                if(vtblA.get(v3,"die")=="die"):

                                                                                                vtblA[v3]=' '.join(vA)

                                                                                else:

                                                                                                vtblA[v3]=vtblA[v3]+" "+' '.join(vA)

                                                                if(e==1):

                                                                                if(vtblB.get(v3,"die")=="die"):

                                                                                                vtblB[v3]=' '.join(vA)

                                                                                else:

                                                                                                vtblB[v3]=vtblB[v3]+" "+' '.join(vA)

                                                                if(e==2):

                                                                                if(vtblC.get(v3,"die")=="die"):

                                                                                                vtblC[v3]=' '.join(vA)

                                                                                else:

                                                                                                vtblC[v3]=vtblC[v3]+" "+' '.join(vA)

'''

 

 

Chapter One. Turn protein fasta into codon optimized nucleotide sequence

 

 

'''

# Get and prepare protein sequences.

vA=[]

for a,b in enumerate(open(vfasta, 'r', encoding="ascii", errors="ignore")):                             #NCBI files often have extraneous BS. Thus the harsh filter.

                b=b.rstrip()

                vA.append(b)

 

# Prepare for EOF.

vA.append(">")                #Ugly Hack.

 

# Transform sequences.

vf1=open(vout, 'w')

for a,b in enumerate(vA):

 

                # Initialize.

                if(">" in b and a==0):

                                v1=b.replace(" ","_")

                                v2=""

 

                if(">" in b and a!=0):

                                vnseq=""

 

                                # Move through last protein sequence (amino acid by amino acid).

                                for c in range(0,len(v2)):

                                                vAA=v2[c:c+1]

 

                                                # Get synonymous codon choices for the present amino acid.

                                                if(vtable_name=="A"):

                                                                v3=vtblA[vAA]

                                                if(vtable_name=="B"):

                                                                v3=vtblB[vAA]

                                                if(vtable_name=="C"):

                                                                v3=vtblC[vAA]

                                                vlist=v3.split(" ")

 

                                                # Randomize synonymous codons and pick one.

                                                random.shuffle(vlist)

                                                vnseq=vnseq+vlist[0]

 

                                if(">" in v1):

                                vf1.write(">"+v1+"|Codon_optimized_for_Ecoli_with_endy_dna20_codons_py_scheme_B_180331"+"\n"+vnseq+"TAATGA"+"\n")

                                else:

                                vf1.write(v1+"|Codon_optimized_for_Ecoli_with_endy_dna20_codons_py_scheme_B_180331"+"\n"+vnseq+"TAATGA"+"\n")

 

                                # Reset some stuff.

                                v1=b.replace(" ","_")

                                v2=""

 

                # Get protein sequence.

                if(">" not in b):

                                v2=v2+b

'''

Should probably write a quality control loop to ensure that no protein accidentally acquires a run of the same 'rare codon or even a run of any 'rare codons.

That would be bad for that proteins expression. But, in reality, it should be density of 'rare codons over some sliding window.

Too much for right now. Have no evidence at my disposal to work from.

'''
