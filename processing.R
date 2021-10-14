library(maftools)
library(tidyverse)

load("data/pancan_nsm_readmaf.RData")

nsm_dat <- nsm.maf@data

summary(as.factor(nsm_dat$Variant_Type))
summary(as.factor(nsm_dat$Variant_Classification))

#only use snp mutations
nsm_snp <- nsm_dat %>%
  filter(Variant_Type == "SNP")

# load prepared data
load(here::here("data/pancan_nsm_readmaf.RData"))
test <- nsm.maf@maf.silent


# read in big unfiltered file
laml <- read.maf(maf = here::here("data/pancan_unfiltered.maf"))

nonsym <- laml@data
silent <- laml@maf.silent

levels(nonsym$Variant_Classification)
levels(nonsym$Variant_Type)
levels(as.factor(nonsym$Consequence))

summary(nonsym$Variant_Classification)
summary


nonsym_noindel <- nonsym %>%
  filter(Variant_Classification %in% c("Nonsense_Mutation"))
levels(nonsym_noindel$Variant_Type)

test <- nonsym_noindel %>%
  filter(Variant_Type %in% c("SNP"))

silent_noindel <- silent %>%
  filter(Variant_Type %in% c("SNP")) %>%
  filter(Variant_Classification %in% "Silent")

levels(as.factor(silent$Variant_Classification))



