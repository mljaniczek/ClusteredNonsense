library(dplyr)
library(data.table)
library(ggplot2)

#This code attempts to categorize the data from the Chang et al paper(https://www.ncbi.nlm.nih.gov/pubmed/26619011) 
#The data will first be filtered to include only nonsense and synonymous mutations, and indels will be excluded
#The genes will then be categorized based on their role in cancer (tumor suppressor, oncogene, other)
#The genes will be further subdivided based on their location on the gene (within 50bp from end of transcript or 50bp upstream of end of transcript)

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
#oncogenes or TSG through a fusion event... For now, I will just exclude those categorized as 
#fusion, as TSG & oncogene, have no category, have three categories (oncogene, TST, fusion)
##Also including Molecular.Genetics categorization of each gene (recessive)
cosmicGenes2 <- cosmicGenes %>% select(Gene.Symbol, Role.in.Cancer, Molecular.Genetics) %>%
  filter(!Role.in.Cancer %in% c("", "fusion")) %>%
  rename("Hugo_Symbol" = "Gene.Symbol") %>%
  mutate(Role.in.Cancer = as.vector(Role.in.Cancer))

#Code that can be used to merge the TSG,fusion and TSG categories, and oncogene,fusion and oncogene categories
cosmicGenes2 <- cosmicGenes2 %>%
  #mutate(Role.in.Cancer = replace(Role.in.Cancer, Role.in.Cancer == "TSG, fusion", "TSG")) %>%
  #mutate(Role.in.Cancer = replace(Role.in.Cancer, Role.in.Cancer == "oncogene, fusion", "oncogene")) %>%
  mutate(Role.in.Cancer = replace(Role.in.Cancer, Role.in.Cancer == "oncogene, TSG, fusion", "oncogene, TSG"))

#Categorizing genes in our data based on role in cancer. 
#Added another column illustrating this  information 

NSM_syn_onco_exon <- left_join(NSM_syn, cosmicGenes2) %>%
  mutate(Role.in.Cancer = replace(Role.in.Cancer, is.na(Role.in.Cancer), "Other"))

##############
##Categorization based on whether the mutation occurs within last exon

#Extracting the exon position and the number of exons per gene and adding to separate columns
NSM_syn_onco_exon$Exon.length <- gsub(".*/", "", NSM_syn_onco_exon$EXON) %>%
  as.numeric
NSM_syn_onco_exon$Exon.pos <- gsub("[-/].*", "", NSM_syn_onco_exon$EXON) %>%
  as.numeric

#Filtering out genes that have only one exon
NSM_syn_onco_exon <- NSM_syn_onco_exon %>%
  filter(Exon.length != 1)

#Adding a column containing boolean values corresponding to whether mutation is within last exon
NSM_syn_onco_exon$In.last.exon <- 
  NSM_syn_onco_exon$Exon.length == NSM_syn_onco_exon$Exon.pos

#saving file. Can later load this file to do any further analysis
save(NSM_syn_onco_exon, file = "NSM_syn_onco_exon.RData")

fisher.test(table(NSM_syn_onco_exon$Variant_Classification, NSM_syn_onco_exon$In.last.exon))

#Fisher exact test analysis to determine enrichment of NSMs in the last exon of genes
#compared to synonymous mutations
TSG <- NSM_syn_onco_exon %>% filter(Role.in.Cancer == "TSG")
fisher.test(table(TSG$Variant_Classification, TSG$In.last.exon))

oncogene <- NSM_syn_onco_exon %>% filter(Role.in.Cancer == "oncogene")
fisher.test(table(oncogene$Variant_Classification, oncogene$In.last.exon))

TSG_dom <- TSG %>% filter(Molecular.Genetics == "Dom")
TSG_Dom_fisher <- fisher.test(table(TSG_dom$Variant_Classification, TSG_dom$In.last.exon)[, c(2,1)])

TSG_rec <- TSG %>% filter(Molecular.Genetics == "Rec")
TSG_rec_fisher <- fisher.test(table(TSG_rec$Variant_Classification, TSG_rec$In.last.exon)[, c(2,1)])

oncogene_dom <- oncogene %>% filter(Molecular.Genetics == "Dom")
oncogene_dom_fisher <- fisher.test(table(oncogene_dom$Variant_Classification, oncogene_dom$In.last.exon)[, c(2,1)])

#Only one oncogene fell under the recessive category
#oncogene_rec <- oncogene %>% filter(Molecular.Genetics == "Rec")
#fisher.test(table(oncogene_rec$Variant_Classification, oncogene_rec$In.last.exon))

TSG_fus <- NSM_syn_onco_exon %>% filter(Role.in.Cancer == "TSG, fusion")
TSG_fus_Dom_fisher <- fisher.test(table(TSG_fus$Variant_Classification, TSG_fus$In.last.exon), c(2,1))

onco_fus <- NSM_syn_onco_exon %>% filter(Role.in.Cancer == "oncogene, fusion")
onco_fus_fisher <- fisher.test(table(onco_fus$Variant_Classification, onco_fus$In.last.exon), c(2,1))

other <- NSM_syn_onco_exon %>%filter(Role.in.Cancer == "Other")
other_fisher <- fisher.test(table(other$Variant_Classification, other$In.last.exon)[, c(2,1)])

OR <- c(TSG_Dom_fisher$estimate, TSG_rec_fisher$estimate, oncogene_dom_fisher$estimate, other_fisher$estimate)
CIlower <- c(TSG_Dom_fisher$conf.int[1], TSG_rec_fisher$conf.int[1], oncogene_dom_fisher$conf.int[1], other_fisher$conf.int[1])
CIupper <- c(TSG_Dom_fisher$conf.int[2], TSG_rec_fisher$conf.int[2], oncogene_dom_fisher$conf.int[2], other_fisher$conf.int[2])
pvalue <- c(TSG_Dom_fisher$p.value, TSG_rec_fisher$p.value, oncogene_dom_fisher$p.value, other_fisher$p.value)
labels <- c("TSG", "TSG", "Oncogene", "Other")
Dominance <- c("Dominant", "Recessive", "Dominant", "")

oddsRatio <- data.frame(labels, OR, CIlower, CIupper, pvalue, Dominance)

ggplot(oddsRatio, aes(x = oddsRatio$label, y = oddsRatio$OR, color = Dominance)) +
  geom_point() + 
  geom_errorbar(aes(ymax = CIupper, ymin = CIlower, width = 0.1)) +
  geom_hline(aes(yintercept = 1), size = .25, linetype = "dashed") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("Odds ratio") +
  xlab("") + 
  ggtitle("Depletion of nonsense mutations in the last exon") + 
  geom_text(aes(label=pvalue2),hjust=-0.1, vjust=0, size = 3, show.legend = FALSE)

TSG_dom_summary <- TSG %>% 
  filter(Variant_Classification == "Nonsense_Mutation", Molecular.Genetics == "Dom") %>%
  group_by(Hugo_Symbol) %>%
  summarise(Total = n(), In.last.exon = sum(In.last.exon == TRUE))
