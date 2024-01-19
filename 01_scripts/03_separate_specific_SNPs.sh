#!/bin/bash
#03_separate_specific_SNPs.sh
#srun -c 1 --mem 10G -p small --time 1-00:00:00 -J 03_separate_SNPs -o 03_separate_SNPs_%j.log ./01_scripts/03_separate_specific_SNPs.sh &

# SNP should be .gz but list without that ending here
SNP="03_SNP_data/all_maf0.05_pctind0.75_maxdepth25.bed"
SNPCOV="03_SNP_data/all_pctind0.75_maxdepth25.bed.gz"
SPECSNPs="05_specific_SNPs"

#SNP file columns are: chr, start, stop, major, minor, reference, probability, num individuals

## all C and G sites in genome (from Claire's unfiltered SNP data)
#gunzip -c "$SNPCOV" |
#awk '($6 == "C" || $6 == "G"){print $1, $2, $3, $4, $5, $6, $7, $8}' - > "$SPECSNPs"/all_pctind0.75_maxdepth25_all_covered_Cs_Gs.bed
#awk -v OFS="\t" 'NR {print $1, $2, $3}' "$SPECSNPs"/all_pctind0.75_maxdepth25_all_covered_Cs_Gs.bed > "$SPECSNPs"/all_pctind0.75_maxdepth25_all_covered_Cs_Gs_3col.bed


gunzip -k "$SNP".gz

### C/T and G/A SNPs
# take SNPS where ref is either C or G, and only print C/T or A/G SNPs
awk '($6 == "C" || $6 == "G") &&
(($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C") || ($4 == "G" && $5 == "A") || ($4 == "A" && $5 == "G")){print $1, $2, $3, $4, $5, $6, $7, $8
    }' "$SNP" > "$SPECSNPs"/all_maf0.05_pctind0.75_maxdepth25_CT_AG_SNPs.bed
awk -v OFS="\t" 'NR {print $1, $2, $3}' "$SPECSNPs"/all_maf0.05_pctind0.75_maxdepth25_CT_AG_SNPs.bed > "$SPECSNPs"/all_maf0.05_pctind0.75_maxdepth25_CT_AG_SNPs_3col.bed


### C/G and G/C SNPs
# take SNPS where ref is either C or G, and only print C/G SNPs (C/G OR G/C)
awk '($6 == "C" || $6 == "G") &&
(($4 == "C" && $5 == "G") || ($4 == "G" && $5 == "C") ){print $1, $2, $3, $4, $5, $6, $7, $8
    }' "$SNP" > "$SPECSNPs"/all_maf0.05_pctind0.75_maxdepth25_CG_SNPs.bed
awk -v OFS="\t" 'NR {print $1, $2, $3}' "$SPECSNPs"/all_maf0.05_pctind0.75_maxdepth25_CG_SNPs.bed > "$SPECSNPs"/all_maf0.05_pctind0.75_maxdepth25_CG_SNPs_3col.bed


### C/A and T/G SNPs
# take SNPS where ref is either C or G, and only print C/A SNPs (C/A OR T/G)
awk '($6 == "C" || $6 == "G") &&
(($4 == "C" && $5 == "A") || ($4 == "A" && $5 == "C") || ($4 == "G" && $5 == "T") || ($4 == "T" && $5 == "G")){
    	print $1, $2, $3, $4, $5, $6, $7, $8
    }' "$SNP" > "$SPECSNPs"/all_maf0.05_pctind0.75_maxdepth25_CA_GT_SNPs.bed
awk -v OFS="\t" 'NR {print $1, $2, $3}' "$SPECSNPs"/all_maf0.05_pctind0.75_maxdepth25_CA_GT_SNPs.bed > "$SPECSNPs"/all_maf0.05_pctind0.75_maxdepth25_CA_GT_SNPs_3col.bed