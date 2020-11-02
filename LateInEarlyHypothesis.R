setwd("/Users/nhansen/Bell_EEC")

# This analysis re-calculates q-values (Benjamini-Hochberg FDRs) for 
# p-values calculated by MutSigCV from 270 early stage tumors for 12 
# genes that had been discovered to be significantly mutated in late-stage 
# tumors.
#
# Nancy F. Hansen
# October 30, 2020

# examine genes significantly mutated in the late group to see if they are SMG in the early group:
earlymutsigresults <- read.table("sigearly_in_late/early_sample_mutsig_results.tdf", header=TRUE)
lateresults <- read.table("sigearly_in_late/late_sample_mutsig_results.tdf", header=TRUE)

latesigresults <- lateresults[lateresults$q.value<=0.1,]
latesiggenes <- latesigresults[,"gene"]
earlysigtests <- earlymutsigresults[!is.na(match(earlymutsigresults$gene, latesiggenes)),] # only the 12 late SMGs
earlysigtests$new12testq.value <- p.adjust(earlysigtests$p.value, method="BH") # PAX6 and KLF3 not SMGs, correcting for 12 tests


