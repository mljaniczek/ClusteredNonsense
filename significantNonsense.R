###########################################################
# Testing for significant nonsense mutations positionally #
#Author: Margie Hannum                                    #
###########################################################

setwd("/Users/TinyDragon/Desktop/Data")
load("NSM_syn_onco.RData") #Amr's subset of nonsense and synonymous mutations, excluding indels, and categorized via Chang et al paper

NSM_onco <- NSM_syn_onco %>% filter(Variant_Classification == "Nonsense_Mutation")
syn_onco <- NSM_syn_onco %>% filter(Variant_Classification == "Silent")

library(MASS) 
library(dplyr)

#Genome-wide quick and dirty fishers-exact test
nsmSynTab <- table(NSM_syn_onco$Variant_Classification, NSM_syn_onco$Within_50bp)
nsmSynTab
fisher.test(nsmSynTab)

#Examining just TP53 last 50 bp T/F
tp53 <- subset(NSM_syn_onco, NSM_syn_onco$Hugo_Symbol == "TP53")
nsmSynTabtp <- table(tp53$Variant_Classification, tp53$Within_50bp)
nsmSynTabtp
fisher.test(nsmSynTab)

#Comparing a certain percent of gene
tp53$Pos_Percent <- 
  (tp53$CDS_length - tp53$CDS_pos)/tp53$CDS_length
nsmSynTabtpPos <- table(tp53$Variant_Classification, (tp53$Pos_Percent >= 0.40))
nsmSynTabtpPos
tp_p <- fisher.test(nsmSynTabtpPos)

#Testing for significant nonsense mutation rate vs background rate of silent mutations


#Getting smaller subset of columns of data so loop can go faster
features = c("Hugo_Symbol", "CDS_length", "CDS_pos", "Variant_Classification", "Within_50bp")
nsmselect <- select(NSM_syn_onco, features)
#Subsetting just genes that have both nonsense and silent mutations
nsmselect <- subset(nsmselect, nsmselect$Hugo_Symbol %in% syn_onco$Hugo_Symbol)
nsmselect <- subset(nsmselect, nsmselect$Hugo_Symbol %in% NSM_onco$Hugo_Symbol)
#Create new vector of the unique remaining genenames
genename <-c(unique(nsmselect$Hugo_Symbol)) #16411 common genes between Nonsense and silent

#Fisher exact test on nonsense vs silent mutations on first half vs second half
start.time = Sys.time()
p = c()
hugo = c()
#for (i in 1108:1109){   #test with a few
for (i in 1:length(genename)){  #for eventual full set
  tryCatch({
    gene = genename[i]
    genesubset <- subset(nsmselect, nsmselect$Hugo_Symbol == gene)
    genesubset$Pos_Percent <- 
    (genesubset$CDS_length - genesubset$CDS_pos)/genesubset$CDS_length
    genesubsetPos <- table(genesubset$Variant_Classification, (genesubset$Pos_Percent >= 0.40))
    genesubsetTest <- fisher.test(genesubsetPos)
    p <- c(p, genesubsetTest$p.value)
    hugo <- c(hugo, gene)},
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}
end.time = Sys.time()
end.time - start.time

#Adjust p-value using Benjamini-Hochberg FDR correction 
p_bh <- data.frame(p.adjust(p, method = "BH")) #where p is vector of p-values
p_bh$gene <- hugo
#Now make subset of p-val under FDR cutoff of 0.1
p_bh_fdr1 <- subset(p_bh, p_bh$p.adjust.p..method....BH.. <= 0.1)

save(p_bh, file = "fdr_fisher_50percent_tot.RData")
save(p_bh_fdr1, file = "fdr_fisher_50percent_select1.RData")

#Results after running on 60/40 ratio
p_bh60 <- data.frame(p.adjust(p, method = "BH")) #where p is vector of p-values
p_bh60$gene <- hugo
#Now make subset of p-val under FDR cutoff of 0.1
p_bh60_fdr1 <- subset(p_bh60, p_bh60$p.adjust.p..method....BH.. <= 0.1)

save(p_bh60, file = "fdr_fisher_60percent_tot.RData")
save(p_bh60_fdr1, file = "fdr_fisher_60percent_select1.RData")

#Results after running on 40/60 ratio
p_bh40 <- data.frame(p.adjust(p, method = "BH")) #where p is vector of p-values
p_bh40$gene <- hugo
#Now make subset of p-val under FDR cutoff of 0.1
p_bh40_fdr1 <- subset(p_bh40, p_bh40$p.adjust.p..method....BH.. <= 0.1)

save(p_bh40, file = "fdr_fisher_40percent_tot.RData")
save(p_bh40_fdr1, file = "fdr_fisher_40percent_select1.RData")


#Loop for last 50 bp (TIME CONSUMING!!!)
start.time = Sys.time()
p = c()
hugo = c()
#for (i in 13849:length(genename)){   #test with a few
for (i in 1:length(genename)){  #for eventual full set
  tryCatch({
    gene = genename[i]
    genesubset <- subset(nsmselect, nsmselect$Hugo_Symbol == gene)
    genesubsetTab <- table(genesubset$Variant_Classification, genesubset$Within_50bp)
    genesubsetTest <- fisher.test(genesubsetTab)
    p <- c(p, genesubsetTest$p.value)
    hugo <- c(hugo, gene)},
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}
end.time = Sys.time()
end.time - start.time

#Adjust p-value using Benjamini-Hochberg FDR correction 
#p_bh_50bp <- data.frame(p.adjust(p, method = "BH")) #where p is vector of p-values
#p_bh_50bp$gene <- hugo

#p_bh_50bp_fdr1 <- subset(p_bh_50bp, p_bh_50bp$p.adjust.p..method....BH.. <= 0.1)

p_bh_50bptot <- data.frame(p.adjust(p, method = "BH")) #where p is vector of p-values
p_bh_50bptot$gene <- hugo
#Now make subset that is under FDR cutoff of 0.1
p_bh_50bpfdr <- subset(p_bh_50bptot, p_bh_50bptot$p.adjust.p..method....BH.. <= 0.1)

save(p_bh_50bpfdr, file = "fdr_50bp_nonsyn.RData")
save(p_bh_50bptot, file = "fdrtot_50bp_nonsyn.RData")



###Examining tumor types per top 8 genes
select <- c("CDKN2A", "EED", "SMAD2", "PIK3R1", "APC", "RB1", "CASP8", "PTEN")

nsm.select <- NSM_onco %>% filter(Hugo_Symbol %in% c("CDKN2A", "EED", "SMAD2", "PIK3R1", "APC", "RB1", "CASP8", "PTEN"))
# 1039 samples across 29 Tumor types have nonsense mutations in our top 8 genes

tumor = c()
count = c()
#for (i in 1:length(select)){
for (i in 1:8){
  tryCatch({
    gene = select[i]
    genesub = subset(nsm.select, nsm.select$Hugo_Symbol == gene)
    count.i = c()
    for (j in 1:length(unique(nsm.select$TUMORTYPE))){
      tumor.j = unique(nsm.select$TUMORTYPE)[j]
      count.i = c(count.i, length(subset(genesub, genesub$TUMORTYPE == tumor.j)$TUMORTYPE))
    }
    count = cbind(count, count.i)
    },
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}

rownames(count) = unique(nsm.select$TUMORTYPE)
colnames(count) = select

##
mat <- matrix(c(789, 4399, 1340, 10084), nrow=2, ncol = 2, byrow = TRUE)
fisher.test(mat)

