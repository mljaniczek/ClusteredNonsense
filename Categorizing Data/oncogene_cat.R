library(dplyr)
library(data.table)
library(ggplot2)

#This code attempts to categorize the data from the Chang et al paper(https://www.ncbi.nlm.nih.gov/pubmed/26619011) 
#The data will first be filtered to include only nonsense and synonymous mutations, and indels will be excluded
#The genes will then be categorized based on their role in cancer (tumor suppressor, oncogene, other)
#The mutations will be further subdivided based on their location on the gene (within 50bp from end of transcript or 50bp upstream of end of transcript)

###I saved the output of the filtering and categorization steps into an RData file that can be loaded direclty

pancan <- fread("pancan_unfiltered.maf")
dim(pancan)

##############

#Filter data to include both nonsense of synonymous mutations (listed under Variant_Classification)
NSM_syn <- pancan %>% filter(Variant_Classification %in% c("Nonsense_Mutation", "Silent"), 
                             !Variant_Type %in% c("DEL", "INS"))

##############

##Categorization based on role in cancer##

#A list of census genes with roles in cancer was downloaded from cosmic(http://cancer.sanger.ac.uk/census)

#Read file containing census genes
cosmicGenes <- read.csv("Census_allThu Nov 16 16_40_30 2017.csv")

##This data contains two "Tiers" of genes. Tier 1 contains established driver genes 
#and Tier 2 genes have less but emerging evidence supporting their role in cancer
#Should we exclude Tier 2 genes?

#Summary of cancer role status for each gene
table(cosmicGenes$Role.in.Cancer)

##Fusion refers to fusion of one gene to another. Typically happens with oncogene, where
#a fusion event results in a protein product that is more activate
##Fusion events occur through translocations, which are mechanistically different from
#point mutations. I can't decide whether to exclude genes that typically become
#oncogenes or TSG through a fusion event... For now I will exclude those that have no
#category, but I kept a code the excludes several other categories commented out

cosmicGenes2 <- cosmicGenes %>% select(Gene.Symbol, Role.in.Cancer) %>%
      filter(Role.in.Cancer != "") %>%
      #filter(!Role.in.Cancer %in% c("", "fusion", "oncogene, TSG, fusion", "oncogene, TSG")) %>%
      rename("Hugo_Symbol" = "Gene.Symbol") %>%
      mutate(Role.in.Cancer = as.vector(Role.in.Cancer))

#Code that can be used to merge the TSG,fusion and TSG categories, and oncogene,fusion and oncogene categories
#cosmicGenes2 <- cosmicGenes2 %>%
    #mutate(Role.in.Cancer = replace(Role.in.Cancer, Role.in.Cancer == "TSG, fusion", "TSG")) %>%
    #mutate(Role.in.Cancer = replace(Role.in.Cancer, Role.in.Cancer == "oncogene, fusion", "oncogene"))
      
  
#Categorizing genes in our data based on role in cancer. 
#Added another column illustrating this  information 

NSM_syn_onco <- left_join(NSM_syn, cosmicGenes2) %>%
  mutate(Role.in.Cancer = replace(Role.in.Cancer, is.na(Role.in.Cancer), "Other"))

##############
##Categorization based on whether the mutation occurs within last 50bp of transcript

#Extracting the position and the length of the coding region and adding to separate columns
NSM_syn_onco$CDS_length <- gsub(".*/", "", NSM_syn_onco$CDS_position) %>%
  as.numeric
NSM_syn_onco$CDS_pos <- gsub("[-/].*", "", NSM_syn_onco$CDS_position) %>%
  as.numeric

#Adding a column containing boolean values corresponding to whether mutation is within 50bp
NSM_syn_onco$Within_50bp <- 
      (NSM_syn_onco$CDS_length - NSM_syn_onco$CDS_pos) <= 50

#saving file. Can later load this file to do any further analysis
save(NSM_syn_onco, file = "NSM_syn_onco.RData")

##############


#Separate NSM and synonymous mutations

NSM_onco <- NSM_syn_onco %>% filter(Variant_Classification == "Nonsense_Mutation")
syn_onco <- NSM_syn_onco %>% filter(Variant_Classification == "Silent")

#Contingency table 
table(NSM_onco$Role.in.Cancer, NSM_onco$Within_50)
