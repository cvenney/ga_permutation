#!/bin/bash
# 00_list_CpGs_in_ref_genome.sh
# list of all CpG sites in reference genome in bed format

GENOME="02_genome/genome.fasta"

# requires fastaRegexFinder
#wget https://github.com/dariober/bioinformatics-cafe/blob/master/fastaRegexFinder/fastaRegexFinder.py?raw=true -O fastaRegexFinder.py
#chmod a+x fastaRegexFinder.py

python fastaRegexFinder.py -f "$GENOME" -r CG --noreverse > 02_genome/all_CpG.bed