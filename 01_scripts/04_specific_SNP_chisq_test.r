#!/usr/bin/env Rscript
#04_specific_SNP_perm_test.R
# conda activate R411 before use
#srun -c 1 --mem 80G -p small --time 1-00:00:00 -J 04_specific_SNPs -o 04_specific_SNPs_%j.log Rscript ./01_scripts/04_specific_SNP_chisq_test.r &

library(regioneR)

allCsGs <- toGRanges("05_specific_SNPs/all_pctind0.75_maxdepth25_all_covered_Cs_Gs_3col.bed")
cg <- toGRanges("04_permutation_files/all_CpG_covered_by_SNPs_nonamecheck_uniques_3col.bed")

ct_snp <- toGRanges("05_specific_SNPs/all_maf0.05_pctind0.75_maxdepth25_CT_AG_SNPs_3col.bed")
cg_snp <- toGRanges("05_specific_SNPs/all_maf0.05_pctind0.75_maxdepth25_CG_SNPs_3col.bed")
ca_snp <- toGRanges("05_specific_SNPs/all_maf0.05_pctind0.75_maxdepth25_CA_GT_SNPs_3col.bed")


### C/T and A/G SNPs
print("####################################### C/T and A/G SNP chi-sq test #######################################")
nCpG <- length(cg) * 2
print(paste("number of CpG site nucleotides covered by SNP data =", nCpG))

ct_nSNP <- length(ct_snp)
print(paste("number of C/T and A/G SNPs =", ct_nSNP))

nOverlaps_ct <- numOverlaps(cg, ct_snp)
print(paste("number of SNP-CpG overlaps =", nOverlaps_ct))

nCsGs <- length(allCsGs)
print(paste("number of Cs and Gs covered by SNP data =", nCsGs))

cont_tab<-rbind(c(nCpG,nOverlaps_ct),
                  c(nCsGs,ct_nSNP))
print(paste("contingency table", cont_tab))

x2 <- chisq.test(cont_tab, correct=FALSE)
print(paste("chi-sq results:"))
print(x2)

print(paste("chi-sq expected values"))
print(x2$expected)


### C/G SNPs
print("####################################### C/G SNP chi-sq test #######################################")

cg_nSNP <- length(cg_snp)
print(paste("number of C/G SNPs =", cg_nSNP))

nOverlaps_cg <- numOverlaps(cg, cg_snp)
print(paste("number of SNP-CpG overlaps =", nOverlaps_cg))

cont_tab<-rbind(c(nCpG,nOverlaps_cg),
                  c(nCsGs,cg_nSNP))
print(paste("contingency table", cont_tab))

x2 <- chisq.test(cont_tab, correct=FALSE)
print(paste("chi-sq results:"))
print(x2)

print(paste("chi-sq expected values"))
print(x2$expected)


### C/A and G/T SNPs
print("####################################### C/A and G/T SNP chi-sq test #######################################")

ca_nSNP <- length(ca_snp)
print(paste("number of C/A SNPs =", ca_nSNP))

nOverlaps_ca <- numOverlaps(cg, ca_snp)
print(paste("number of SNP-CpG overlaps =", nOverlaps_ca))

cont_tab<-rbind(c(nCpG,nOverlaps_ca),
                  c(nCsGs,ca_nSNP))
print(paste("contingency table", cont_tab))

x2 <- chisq.test(cont_tab, correct=FALSE)
print(paste("chi-sq results:"))
print(x2)

print(paste("chi-sq expected values"))
print(x2$expected)
