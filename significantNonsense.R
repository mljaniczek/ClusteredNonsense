###########################################################
# Testing for significant nonsense mutations positionally #
#Author: Margie Hannum                                    #
###########################################################

setwd("/Users/TinyDragon/Desktop/Data")
load("NSM_syn_onco.RData") #Amr's subset of 

NSM_onco <- NSM_syn_onco %>% filter(Variant_Classification == "Nonsense_Mutation")
syn_onco <- NSM_syn_onco %>% filter(Variant_Classification == "Silent")

library(MASS) 

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
nsmSynTabtpPos <- table(tp53$Variant_Classification, (tp53$Pos_Percent >= 0.60))
nsmSynTabtpPos
tp_p <- fisher.test(nsmSynTabtpPos)

#Getting smaller subset of columns of data so loop can go faster
features = c("Hugo_Symbol", "CDS_length", "CDS_pos", "Variant_Classification", "Within_50bp")
nsmselect <- select(NSM_syn_onco, features)
nsmselect <- subset(nsmselect, nsmselect$Hugo_Symbol %in% syn_onco$Hugo_Symbol)
nsmselect <- subset(nsmselect, nsmselect$Hugo_Symbol %in% NSM_onco$Hugo_Symbol)
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
    genesubsetPos <- table(genesubset$Variant_Classification, (genesubset$Pos_Percent >= 0.50))
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
p_bh_fdr1 <- subset(p_bh, p_bh$p.adjust.p..method....BH.. <= 0.1)

save(p_bh, file = "fdr_fisher_50percent_tot.RData")
save(p_bh_fdr1, file = "fdr_fisher_50percent_select1.RData")


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

p_bh_50bpfdr <- subset(p_bh_50bptot, p_bh_50bptot$p.adjust.p..method....BH.. <= 0.1)

save(p_bh_50bpfdr, file = "fdr_50bp_nonsyn.RData")
save(p_bh_50bptot, file = "fdrtot_50bp_nonsyn.RData")

