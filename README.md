# FreeGenes

This is the FreeGenes repository, which holds code that links the online submission forms to the FreeGenes Bionet. 

# Pipeline

## retrieve.py
"retrieve.py" is a script that takes information from an online spreadsheet and the Bionet with entries from it.

This program first downloads and creates a pandas dataframe based off of a Google Sheet. This Google Sheet can be created by a variety of ways, including a Google Form or by direct entry. We use a Wufoo form combined with an Automate.io bot to automatically fill the Google Sheet. Depending on which method you use, the way "Genbank File" and "CSV File" (for bulk uploads) are downloaded has to be changed. 

There are 2 different upload methods: Single submission and Bulk submission. Though they need to be handled in slightly different ways, both are written as a FreeGene obkect

After retrieval of sequences, CDSs are optimized to remove restriction enzyme sites (as defined in "optimize.py"). If a non-CDS has BbsI, BsaI, or BtgZI, it is rejected and an email is sent to the contributor (IMPLEMENT). All this data is loaded into a FreeGene object for each part submitted. 

If all checks(IMPLEMENT) are completed, then the FreeGene object is written to a Bionet node. From there, genes can be loaded into the OpenFoundry for processing into physical molecules.

## config.yaml

## template.json

## other scripts

 
