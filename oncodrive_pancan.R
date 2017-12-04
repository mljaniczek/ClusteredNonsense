library(dplyr)
library(data.table)
library(ClusteredMutations)
library(ggplot2)

#install and load maftools
source("https://bioconductor.org/biocLite.R")
biocLite("maftools")
library(maftools)

#oncodrive with our pancan data
#Filter for nonsense mutations but keep all columns 

#######################################
# Oncodrive w with Amr's filtered data #
#######################################
setwd("/Users/TinyDragon/Desktop/Data")
load("NSM_syn_onco.RData") #Amr's subset of nonsense and synonymous mutations, excluding indels, and categorized via Chang et al paper

##Clustering using all gene types
#Use read.maf from maftools package to prepare MAF class object
nsm.maf <- read.maf(maf = NSM_syn_onco, useAll = TRUE)
#saved maf to new file for later use
save(nsm.maf, file="pancan_nsm_new_readmaf.RData")
load("pancan_nsm_new_readmaf.RData")
#Run oncodrive function on maf object (took 20 minutes)
#With bgEstimate = FALSE, indicate that function should use background mutation rate from COSMIC
nsm.sig15 <- oncodrive(maf = nsm.maf, AACol = 'Amino_Acid_Change', 
                     minMut = 15, pvalMethod = 'zscore', bgEstimate = TRUE)
#estimated bg scores for synonymous variants mean 0.23138, estimated bg sd: 0.12427

#Select results from oncodriveClust that have >=5 mutations in the clusters
nsm.sig.sel <- subset(nsm.sig, nsm.sig$muts_in_clusters >= 5) 
nsm.sig.sel15 <- subset(nsm.sig15, nsm.sig15$muts_in_clusters >= 5) 

#Save result for later use
saveRDS(nsm.sig.sel, file = "pancan_nsm_oncodrive_sig_sel.rds")
load("pancan_nsm_oncodrive_sig_sel.rds") #load if necessary

#basic plut function from MAF
plotOncodrive(res = nsm.sig.sel15, fdrCutOff = 0.05, useFraction = TRUE, labelSize = 3) 

#Subset oncodriveCLUST results further 
nsm.select <- subset(nsm.sig.sel, nsm.sig.sel$clusterScores > 0.5)
nsm.select <- subset(nsm.sig.sel, nsm.sig.sel$fdr < 0.05)
nsm.select15 <- subset(nsm.sig.sel15, nsm.sig.sel15$fdr < 0.05)
save(nsm.select15, file="nsm.select15.RData")
#Read file containing census genes, to classify
cosmicGenes <- read.csv("Census_allMon Nov 20 16_36_34 2017.csv")
    
cosmicSelect <- subset(cosmicGenes, cosmicGenes$Gene.Symbol %in% nsm.select15$Hugo_Symbol)

s<- cosmicSelect$Gene.Symbol

role <- c(rep("Other", 37))
tsgind <- c(1, 5, 7, 8, 10, 15, 18, 19, 20, 26, 34, 36)
role <- replace(role, tsgind, "TSG")
role <- replace(role, 24, "Oncogene")

nsm.select15$Role.in.cancer <- as.factor(role)

#Pulled code for plotOncodrive since I wanted to make some modifications
res=nsm.select15
fdrCutOff = 0.05
labelSize=5
res$label = paste(res$Hugo_Symbol, '[',res$clusters,']', sep='')
res$significant = ifelse(test = res$fdr < fdrCutOff, yes = 'sig', no = 'nonsig')
colCode = c('sig' = 'red', 'nonsig' = 'royalblue')
tsgplot = ggplot(data = res, aes(x = fract_muts_in_clusters, y = -log10(fdr), size = clusters, color = significant))+
  geom_point(alpha = 0.9)+cowplot::theme_cowplot(line_size = 1)+theme(legend.position = 'NONE')+scale_color_manual(values = colCode)+
  ggrepel::geom_text_repel(data = res[fdr < fdrCutOff], aes(x = fract_muts_in_clusters, y = -log10(fdr), label = label, size = labelSize), color = 'black')+
  xlab('Fraction of mutations in clusters')+cowplot::background_grid(major = 'xy') + ggtitle('Identification of genes containing nonsense mutation clusters')

#Changing color to class of gene
nsmplot = ggplot(data = res, aes(x = fract_muts_in_clusters, y = -log10(fdr), size = clusters, color = Role.in.cancer))+
  geom_point(alpha = 0.9)+cowplot::theme_cowplot(line_size = 1)+#theme(legend.position = 'NONE')+#scale_color_manual(values = colCode)+
  ggrepel::geom_text_repel(data = res[fdr < fdrCutOff], aes(x = fract_muts_in_clusters, y = -log10(fdr), label = label, size = labelSize), color = 'black')+
  xlab('Fraction of mutations in clusters')+cowplot::background_grid(major = 'xy') + ggtitle('Identification of genes containing nonsense mutation clusters')


summary(nsm.sig.sel)





##Old code
################################################################
# Reading and filtering pancan data - you can just skip this   # 
# and load prepared files below this section!!                 #
################################################################
setwd("/Users/TinyDragon/Desktop/Data")
pancan <- fread("pancan_unfiltered.maf")
nsmAll <- select(pancan, colnames(pancan)) %>%
  filter(Variant_Classification == "Nonsense_Mutation")

nsmAll <- select(pancan, colnames(pancan)) %>%
  filter(Variant_Classification == "Nonsense_Mutation")

nsmsilent <- select(pancan, colnames(pancan)) %>%
  filter(Variant_Classification == "Silent") 

#Use read.maf from maftools package to prepare MAF class object
nsm.maf <- read.maf(maf = nsmAll, useAll = TRUE)
#saved maf to new file for later use
save(nsm.maf, file="pancan_nsm_readmaf.RData")
#Run oncodrive function on maf object (took 20 minutes)
#With bgEstimate = FALSE, indicate that function should use background mutation rate from COSMIC
nsm.sig <- oncodrive(maf = nsm.maf, AACol = 'Amino_Acid_Change', 
                     minMut = 5, pvalMethod = 'zscore', bgEstimate = FALSE)
save(nsm.sig, file = "pancan_nsmsig_oncodrive.RData") #Save result
######################################################################


#########################################################################
# You can  load below prepared files and run the plots if you want! #
#########################################################################

#First set your working directory to wherever you save the files
setwd("<INSERT LOCATION HERE>")
load("pancan_nsm_readmaf.RData") #object nsm.maf which is just the nonsense mutations but after going through the read.maf function to make a maf object
load("pancan_nsmsig_oncodrive.RData") #For nsm.sig - the result of OncodriveClust algorithm

#Plot result of oncodrive using FDR = 0.05,
plotOncodrive(res = nsm.sig, fdrCutOff = 0.05, useFraction = FALSE)
plotOncodrive(res = nsm.sig, fdrCutOff = 0.05, useFraction = TRUE)

#All columns necessary for maftools functions 

#Check summary of oncodrive cluster results across 16502 genes
summary(nsm.sig$Nonsense_Mutation) #Number of nonsense mutations per gene Min 5, Q1 7, Med 10, Mean 14.69, Q3 16, Max 572
summary(nsm.sig$clusters) #Number of clusters per gene Min 1, Q1 1, Med 1, Mean 1.7, Q3 2, Max 45
summary(nsm.sig$muts_in_clusters) #Number of mutations per cluster Min 2, Q1 2, Med 3, Mean 4.63, Q3 5, Max 540
summary(nsm.sig$fract_muts_in_clusters) #Fraction of total mutations in clusters per gene Min 0.04, Q1 0.22, Med 0.31, Mean 0.33, Q3 0.4, Max 1
summary(nsm.sig$fdr) #FDR per gene Min 7.7e-06, Q1 0.55, Med 0.82, Mean 0.688, Q3 0.89, Max 0.966
length(nsm.sig$fdr >= 0.05) ##3694 Genes had FDR >= 0.05


## Using Amr's filtered data, subsetting based on class
tsg <- subset(NSM_syn_onco, NSM_syn_onco$Role.in.Cancer == "TSG") #163 genes, 4760 samples
#Use read.maf from maftools package to prepare MAF class object
tsg.maf <- read.maf(maf = tsg, useAll = TRUE)
#saved maf to new file for later use
save(tsg.maf, file="pancan_nsm_tsg_readmaf.RData")
load("pancan_nsm_tsg_readmaf.RData")
#Run oncodrive function on maf object (took 20 minutes)
#With bgEstimate = FALSE, indicate that function should use background mutation rate from COSMIC
tsg.sig <- oncodrive(maf = tsg.maf, AACol = 'Amino_Acid_Change', 
                     minMut = 5, pvalMethod = 'zscore', bgEstimate = TRUE)
tsg.sig.sel <- subset(tsg.sig, tsg.sig$muts_in_clusters >= 4)
plotOncodrive(res = tsg.sig.sel, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 3) 

res=tsg.sig
fdrCutOff = 0.1
labelSize=3
res$label = paste(res$Hugo_Symbol, '[',res$clusters,']', sep='')
res$significant = ifelse(test = res$fdr < fdrCutOff, yes = 'sig', no = 'nonsig')
colCode = c('sig' = 'red', 'nonsig' = 'royalblue')
tsgplot = ggplot(data = res, aes(x = fract_muts_in_clusters, y = -log10(fdr), size = clusters, color = significant))+
  geom_point(alpha = 0.9)+cowplot::theme_cowplot(line_size = 1)+theme(legend.position = 'NONE')+scale_color_manual(values = colCode)+
  ggrepel::geom_text_repel(data = res[fdr < fdrCutOff], aes(x = fract_muts_in_clusters, y = -log10(fdr), label = label, size = labelSize), color = 'black')+
  xlab('Fraction of mutations in clusters')+cowplot::background_grid(major = 'xy') + ggtitle('TSG Genes with Nonsense Mutation Clusters')

## NOW on TSG Fustion genes
tsg.fus <- subset(NSM_syn_onco, NSM_syn_onco$Role.in.Cancer == "TSG, fusion") #60 genes, 1894 samples
#Use read.maf from maftools package to prepare MAF class object
tsgfus.maf <- read.maf(maf = tsg.fus, useAll = TRUE)
#Run oncodrive function on maf object
#With bgEstimate = FALSE, indicate that function should use background mutation rate from COSMIC
tsgfus.sig <- oncodrive(maf = tsgfus.maf, AACol = 'Amino_Acid_Change', 
                        minMut = 5, pvalMethod = 'zscore', bgEstimate = TRUE)
tsgfus.sig.sel <- subset(tsgfus.sig, tsgfus.sig$muts_in_clusters >=4)
plotOncodrive(res = tsgfus.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 3)

res=tsgfus.sig
fdrCutOff = 0.1
labelSize=3
res$label = paste(res$Hugo_Symbol, '[',res$clusters,']', sep='')
res$significant = ifelse(test = res$fdr < fdrCutOff, yes = 'sig', no = 'nonsig')
colCode = c('sig' = 'red', 'nonsig' = 'royalblue')
tsgfusplot = ggplot(data = res, aes(x = fract_muts_in_clusters, y = -log10(fdr), size = clusters, color = significant))+
  geom_point(alpha = 0.9)+cowplot::theme_cowplot(line_size = 1)+theme(legend.position = 'NONE')+scale_color_manual(values = colCode)+
  ggrepel::geom_text_repel(data = res[fdr < fdrCutOff], aes(x = fract_muts_in_clusters, y = -log10(fdr), label = label, size = labelSize), color = 'black')+
  xlab('Fraction of mutations in clusters')+cowplot::background_grid(major = 'xy') + ggtitle('TSG Fusion Genes with Nonsense Mutation Clusters')


## NOW on Oncogenes
onco <- subset(NSM_syn_onco, NSM_syn_onco$Role.in.Cancer == "oncogene") #101 genes, 2989 samples
#Use read.maf from maftools package to prepare MAF class object
onco.maf <- read.maf(maf = onco, useAll = TRUE)
#Run oncodrive function on maf object
#With bgEstimate = FALSE, indicate that function should use background mutation rate from COSMIC
onco.sig <- oncodrive(maf = onco.maf, AACol = 'Amino_Acid_Change', 
                      minMut = 5, pvalMethod = 'zscore', bgEstimate = TRUE)
onco.sig.sel <- subset(tsgfus.sig, tsgfus.sig$muts_in_clusters >=4)
plotOncodrive(res = onco.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 3)

res=onco.sig
fdrCutOff = 0.1
labelSize=3
res$label = paste(res$Hugo_Symbol, '[',res$clusters,']', sep='')
res$significant = ifelse(test = res$fdr < fdrCutOff, yes = 'sig', no = 'nonsig')
colCode = c('sig' = 'red', 'nonsig' = 'royalblue')
oncoplot = ggplot(data = res, aes(x = fract_muts_in_clusters, y = -log10(fdr), size = clusters, color = significant))+
  geom_point(alpha = 0.9)+cowplot::theme_cowplot(line_size = 1)+theme(legend.position = 'NONE')+scale_color_manual(values = colCode)+
  ggrepel::geom_text_repel(data = res[fdr < fdrCutOff], aes(x = fract_muts_in_clusters, y = -log10(fdr), label = label, size = labelSize), color = 'black')+
  xlab('Fraction of mutations in clusters')+cowplot::background_grid(major = 'xy') + ggtitle('Oncogenes with Nonsense Mutation Clusters')

## NOW on Oncogene fusion
oncofus <- subset(NSM_syn_onco, NSM_syn_onco$Role.in.Cancer == "oncogene, fusion") #136 genes, 2602 samples
#Use read.maf from maftools package to prepare MAF class object
oncofus.maf <- read.maf(maf = oncofus, useAll = TRUE)
#Run oncodrive function on maf object
#With bgEstimate = FALSE, indicate that function should use background mutation rate from COSMIC
oncofus.sig <- oncodrive(maf = oncofus.maf, AACol = 'Amino_Acid_Change', 
                         minMut = 5, pvalMethod = 'zscore', bgEstimate = TRUE)
oncofus.sig.sel <- subset(tsgfus.sig, tsgfus.sig$muts_in_clusters >=4)
plotOncodrive(res = oncofus.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 3)

res=oncofus.sig
fdrCutOff = 0.1
labelSize=3
res$label = paste(res$Hugo_Symbol, '[',res$clusters,']', sep='')
res$significant = ifelse(test = res$fdr < fdrCutOff, yes = 'sig', no = 'nonsig')
colCode = c('sig' = 'red', 'nonsig' = 'royalblue')
oncofusplot = ggplot(data = res, aes(x = fract_muts_in_clusters, y = -log10(fdr), size = clusters, color = significant))+
  geom_point(alpha = 0.9)+cowplot::theme_cowplot(line_size = 1)+theme(legend.position = 'NONE')+scale_color_manual(values = colCode)+
  ggrepel::geom_text_repel(data = res[fdr < fdrCutOff], aes(x = fract_muts_in_clusters, y = -log10(fdr), label = label, size = labelSize), color = 'black')+
  xlab('Fraction of mutations in clusters')+cowplot::background_grid(major = 'xy') + ggtitle('Oncogenes, fusion with Nonsense Mutation Clusters')


## NOW on Other gene
other <- subset(NSM_syn_onco, NSM_syn_onco$Role.in.Cancer == "Other") #18795 genes, 10602 samples
#Use read.maf from maftools package to prepare MAF class object
other.maf <- read.maf(maf = other, useAll = TRUE)
#Run oncodrive function on maf object
#With bgEstimate = FALSE, indicate that function should use background mutation rate from COSMIC
other.sig <- oncodrive(maf = other.maf, AACol = 'Amino_Acid_Change', 
                       minMut = 5, pvalMethod = 'zscore', bgEstimate = TRUE)
other.sig.sel <- subset(tsgfus.sig, tsgfus.sig$muts_in_clusters >=4)
plotOncodrive(res = other.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 3)

res=other.sig
fdrCutOff = 0.1
labelSize=3
res$label = paste(res$Hugo_Symbol, '[',res$clusters,']', sep='')
res$significant = ifelse(test = res$fdr < fdrCutOff, yes = 'sig', no = 'nonsig')
colCode = c('sig' = 'red', 'nonsig' = 'royalblue')
oncofusplot = ggplot(data = res, aes(x = fract_muts_in_clusters, y = -log10(fdr), size = clusters, color = significant))+
  geom_point(alpha = 0.9)+cowplot::theme_cowplot(line_size = 1)+theme(legend.position = 'NONE')+scale_color_manual(values = colCode)+
  ggrepel::geom_text_repel(data = res[fdr < fdrCutOff], aes(x = fract_muts_in_clusters, y = -log10(fdr), label = label, size = labelSize), color = 'black')+
  xlab('Fraction of mutations in clusters')+cowplot::background_grid(major = 'xy') + ggtitle('"Other" genes with Nonsense Mutation Clusters')


