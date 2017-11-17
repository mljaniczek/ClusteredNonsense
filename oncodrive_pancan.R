library(dplyr)
library(data.table)
library(ClusteredMutations)

#install and load maftools
source("https://bioconductor.org/biocLite.R")
biocLite("maftools")
library(maftools)

#Oncodrive with example data from maftools package, if you want to test it
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = laml.maf, useAll = FALSE)
laml.sig <- oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1)

#oncodrive with our pancan data
#Filter for nonsense mutations but keep all columns 
#All columns necessary for maftools functions 

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

#Check summary of oncodrive cluster results across 16502 genes
summary(nsm.sig$Nonsense_Mutation) #Number of nonsense mutations per gene Min 5, Q1 7, Med 10, Mean 14.69, Q3 16, Max 572
summary(nsm.sig$clusters) #Number of clusters per gene Min 1, Q1 1, Med 1, Mean 1.7, Q3 2, Max 45
summary(nsm.sig$muts_in_clusters) #Number of mutations per cluster Min 2, Q1 2, Med 3, Mean 4.63, Q3 5, Max 540
summary(nsm.sig$fract_muts_in_clusters) #Fraction of total mutations in clusters per gene Min 0.04, Q1 0.22, Med 0.31, Mean 0.33, Q3 0.4, Max 1
summary(nsm.sig$fdr) #FDR per gene Min 7.7e-06, Q1 0.55, Med 0.82, Mean 0.688, Q3 0.89, Max 0.966
length(nsm.sig$fdr >= 0.05) ##3694 Genes had FDR >= 0.05


oncoplot(maf = nsm.maf, top = 10, fontSize = 12)

apc.lpop = lollipopPlot(maf = nsm.maf, gene = 'APC', AACol = 'Amino_Acid_Change', showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE)
