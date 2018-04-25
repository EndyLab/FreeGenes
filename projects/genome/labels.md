# Potential labels
{'/spliced_product=', '/db_xref=', '/feature_type=', '/protein_id=', '/EC_number=', '/GenBank_acc=', '/spliced_annotation=', '/translation=', '/strand=', '/Source=', '/codon_start=', '/label=', '/ribosomal_slippage/gene=', '/experiment=', '/gene_synonym=', '/locus_tag=', '/function=', '/transl_table=', '/gene=', '/product=', '/note='}
(ecoli)


#Ecoli specific
- /label refers to the gene name
- /spliced_product specifically refers to RF-2, which uses natural UGA slippage
- /experiment is a PMID link
- /function is the function of the gene
- /spliced annotation
- /ribosomal slip


## Metadata
- /locus_tag is the locus tag 
- /translation
## Database stuff
- /protein_id links to ncbi ("https://www.ncbi.nlm.nih.gov/protein/NP_417760.1")
- /EC_number links to enzyme number ("https://enzyme.expasy.org/EC/6.1.1.2")
- /GenBank_acc links to source ("https://www.ncbi.nlm.nih.gov/nuccore/NC_000913")
## Gene data
- /transl_table
- /Source
- /gene_synonyms
- /gene
- /product
- /note


## Processed separately
- /db_xref are database links (multiple links)


- Use locus tag to http://ecoliwiki.net/colipedia/index.php/b3056




# NOT USEFUL
- /feature_type "CDS"
- spliced_annotation
- strand
- codon_start
- ribosomal_slippage

