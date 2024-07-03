###################################################
#### Script Name: 00_Impute_protein_npx
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to impute missing NPX values for proteomics in UKB 
#### Date : 23/03/2023

rm(list = ls())

## load packages 
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(RNOmni)
# library(glue)
library(DBI)
library(dbplyr)
library(magrittr)
# library(broom)
library(survival)
library(data.table)
library(glue)
library(parallel)
library(caret)
library(glmnet)
library(ROSE)
# library(odbc)

## define panel for imputation 
args <- commandArgs(trailingOnly = T)

panel.in <- as.character(args[1])

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/00_protein_imputation/bin/") 

## get protein data and metadata
proteomics <-  fread(cmd = "zcat /home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072.txt.gz", sep = "\t", data.table = F)

ol.npx <- proteomics %>% 
  filter(ins_index == 0&covid_imaging==F) %>% 
  select(eid_20361, sample_id,assay, olink_id,npx , randomized_baseline) %>% 
  # filter(randomized_baseline==T) %>% 
  mutate(mrc_olink.id=gsub("-","_",paste0(assay,"_",gsub("OID","", olink_id)))) %>% 
  pivot_wider(names_from=mrc_olink.id, id_cols = c(eid_20361, randomized_baseline),values_from=npx) 

protein_list <- proteomics %>% 
  mutate(mrc_olink.id=gsub("-","_",paste0(assay,"_",gsub("OID","", olink_id))))%>% 
  select(mrc_olink.id,panel ) %>% 
  distinct()

protein_list$mrc_olink.id <- gsub("-","_", protein_list$mrc_olink.id)

## read in data 
cov <- as.data.frame(fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", sep = "\t"))

## merge - note this is dropping one participant with no covariates
pheno <- merge(cov[,c(1,11,17)], ol.npx, by="eid_20361")

## rename to olds names to match the rest fo the script
colnames(pheno)[c(2)] <- c("ge_age")

## exclude participants with missing age, sex
pheno <- pheno[complete.cases(pheno[,c("ge_sex", "ge_age")]),]

## exclude proteins with more 40 % missingness
# colSums(is.na(pheno[,protein_list$mrc_olink.id]))[which(colSums(is.na(pheno[,protein_list$mrc_olink.id]))>(51810*.3))]
## none proteins need to be excluded based on missing values 

## exclude individuals with more than 60 % of missingness
# ex.ids <- pheno$eid_20361[which(rowSums(is.na(pheno[,protein_list$mrc_olink.id]))>(1472*0.6))]
# pheno <- pheno %>% 
#   filter(!eid_20361%in%ex.ids)

library(missForest)
# library(doMC)
doParallel::registerDoParallel(cores = 96)

## only use age and sex to impute other missing proteins npx
imp.data <- pheno[,c("ge_age", "ge_sex", protein_list$mrc_olink.id[which(protein_list$panel==panel.in)])]
row.names(imp.data) <- pheno$eid_20361
imp.data$ge_sex <- as.factor(imp.data$ge_sex)

m.forest <- missForest(imp.data, verbose = T, replace = F, parallelize = "forests",
                       # variablewise = T, 
                       ntree = 50)
save(m.forest, file = paste0("../data_input/Imputed_NPX_missForest_",panel.in,".RData"))
m.forest$OOBerror
