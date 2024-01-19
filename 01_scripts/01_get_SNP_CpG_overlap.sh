#!/bin/bash
# 01_get_SNP_CpG_overlap.sh
# srun -c 1 --mem 35G -p small --time 1-00:00:00 -J 01_overlap -o 01_overlap_%j.log ./01_scripts/01_get_SNP_CpG_overlap.sh &

ALLCpG="02_genome/all_CpG.bed"
SNPCOV="03_SNP_data/all_pctind0.75_maxdepth25.mafs.gz"
SNPs="03_SNP_data/all_maf0.05_pctind0.75_maxdepth25.mafs.gz"

module load bedtools

base="$(ls -1 "$SNPCOV" | cut -d "." -f1,2)"
gunzip -c "$base".mafs.gz |
awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1, $2-1, $2, $3, $4, $5, $6, $7}' - > "$base".bed
bedtools intersect -a "$base".bed -b "$ALLCpG" -wb -nonamecheck |
awk '!seen[$12]++' - > 04_permutation_files/all_CpG_covered_by_SNPs_nonamecheck_uniques_with_info.bed
awk -v FS="\t" -v OFS="\t" 'NR {print $9, $10, $11}' 04_permutation_files/all_CpG_covered_by_SNPs_nonamecheck_uniques_with_info.bed > 04_permutation_files/all_CpG_covered_by_SNPs_nonamecheck_uniques_3col.bed

base2="$(ls -1 "$SNPs" | cut -d "." -f1,2,3)"
echo "$base2"
gunzip -c "$base2".mafs.gz |
awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1, $2-1, $2, $3, $4, $5, $6, $7}' - > "$base2".bed
gzip "$base2".bed