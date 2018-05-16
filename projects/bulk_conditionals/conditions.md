# Conditions for FreeGene building
In many cases, a DNA part will require certain conditions in order for it to be built. These can only be added to a bulk upload sheet by adding addition columns with specific titles and conditions. There are 7 sections that must be added: "conditions", "vector", "strain", "selection", "media", "temperature", and "notes".

You will be notified if there are any conditionals to be met when building. 


# Adding conditionals
{"vector": "", "strain": "","selection": [""], "media": "", "temperature": "", "mastermix_notes": ""}

## vector
- Options: List of potential vectors. One will be choosen if many are given. 
- Example: ["Entry2"]
- Acceptable cases ["Entry2"]
This is where you can put a desired vector. The default is whatever default entry vector we are using at the time. Entry1 has been depreciated. At the current time, Entry2 type vectors are the only ones available. 

## strain
- Options: Desired cloning strain. 
- Example: "Top10"
- Acceptable cases ["Top10", "ccdB"]
Strain information. Especially important for ccdB

## selection
- Options: List of antibiotics to add to media
- Example: ["Amp"]
- Acceptable cases ["Amp","Kan","Cam","Tet","Spc","Gnt"]
Antibiotics to add to media. Important when creating vectors. 

## media
- Options: Media to plate on
- Example: "IPTG"
- Acceptable cases - CASE BY CASE BASIS.
Should very rarely be filled. Each case check manually to make proper media.

## temperature
- Options: Temperature to grow at
- Example: "37c"
- Acceptable cases - ["RT","30c","37c"]
For temperature sensitive vectors.

## mastermix notes
- Options: The mastermix formulation
- Example: IN PROGRESS
- Acceptable cases: IN PROGRESS
Eventually, numbered tubes of mastermix are desired. 




