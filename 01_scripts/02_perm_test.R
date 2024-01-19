#!/usr/bin/env Rscript
# regioneR to count SNP-CpG overlaps, permutation test, using conda env R11
# need to output results to logfile!!
#srun -c 1 --mem 10G -p small --time 1-00:00:00 -J 02_permtest -o 02_permtest_%j.log Rscript ./01_scripts/02_perm_test.R &

library(regioneR)

#use "gunzip -c 03_SNP_data/all_pctind0.75_maxdepth25.bed.gz | wc -l" or "wc -l 03_SNP_data/all_pctind0.75_maxdepth25.bed" for genome length
GEN_LENGTH=1998994058
NUM_PERM=10000
cg <- toGRanges("04_permutation_files/all_CpG_covered_by_SNPs_nonamecheck_uniques_3col.bed")
snp <- toGRanges("03_SNP_data/all_maf0.05_pctind0.75_maxdepth25.bed.gz")

NUM_CpGs <- length(cg)*2
NUM_SNPs <- length(snp)
SNPs_IN_CpGs <- suppressWarnings(numOverlaps(cg, snp))
PROP_GEN_CpGs <- ((NUM_CpGs) / GEN_LENGTH)

print(paste("number of CpG site nucleotides covered by SNP data =", NUM_CpGs))

print(paste("number of SNPs = ", NUM_SNPs))

print(paste("number of SNPs in CpGs", SNPs_IN_CpGs))

print(paste("proportion of genome made up of CpGs", PROP_GEN_CpGs))


cont_tab<-rbind(c(NUM_SNPs,SNPs_IN_CpGs),
                  c(GEN_LENGTH,NUM_CpGs))
print(paste("contingency table", cont_tab))

x2 <- chisq.test(cont_tab, correct=FALSE)
print(paste("chi-sq results:"))
print(x2)


# permutation test
# count (NUM_PERM) times how many SNPs fell randomly within a region of the genome that was PROP_GEN_CpGs
# then sum the number of times that this random number of SNPs was more than SNPs_IN_CpGs
#x = NULL; for (i in 1:NUM_PERM){print(i);x=c(x, sum(runif(NUM_SNPs) < PROP_GEN_CpGs))}; sum(x > SNPs_IN_CpGs)
#print(" / number of permutations = pval  ")
