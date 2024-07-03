###################################################
#### Script Name: 01_save_exclusions_by_missing_rate
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script determine individual exclusion based on missing rate of NPX values
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
dim(ol.npx)

protein_list <- proteomics %>% 
  mutate(mrc_olink.id=gsub("-","_",paste0(assay,"_",gsub("OID","", olink_id))))%>% 
  select(mrc_olink.id,panel , UniProt) %>% 
  distinct()

protein_list$mrc_olink.id <- gsub("-","_", protein_list$mrc_olink.id)
write.table(protein_list, file = "../data_input/protein_labels_by_panel.txt", sep = "\t", row.names = F)

## read in data 
cov <- as.data.frame(fread("/home/uk_biobank/proteomics/data/processed/ukb_ppp_1536_and_3072_covariate_annotations.txt", sep = "\t"))
## merge - note this is dropping one participant with no covariates
pheno <- merge(cov[,c(1,11,17)], ol.npx, by="eid_20361")
## rename to olds names to match the rest fo the script
colnames(pheno)[c(2)] <- c("ge_age")

## exclude participants with missing age, sex
pheno <- pheno[complete.cases(pheno[,c("ge_sex", "ge_age")]),]
dim(pheno)
write.table(pheno[,c("eid_20361","ge_sex", "ge_age", "randomized_baseline")], file = "../data_input/temp_imputation_ids.txt", sep = "\t", row.names = F)

## none proteins need to be excluded based on missing values
ex.mis.rate <- list()
ex.mis.rate$proteins <- names(which(colSums(is.na(pheno[,protein_list$mrc_olink.id]))>(nrow(pheno)*.25)))

## save missing rate
ex.mis.rate$protein.mis.rate <- (colSums(is.na(pheno[,protein_list$mrc_olink.id]))[which(colSums(is.na(pheno[,protein_list$mrc_olink.id]))>(nrow(pheno)*.25))]*100)/nrow(pheno)

## individual's missing rate
ex.ids <- pheno$eid_20361[which(rowSums(is.na(pheno[,protein_list$mrc_olink.id]))>(length(protein_list$mrc_olink.id)*0.25))]
length(ex.ids)
ex.mis.rate$individuals <- ex.ids

ex.mis.rate$individuals.mis.rate <- (rowSums(is.na(pheno[,protein_list$mrc_olink.id]))[which(rowSums(is.na(pheno[,protein_list$mrc_olink.id]))>(length(protein_list$mrc_olink.id)*0.25))]*100)/length(protein_list$mrc_olink.id)
save(ex.mis.rate, file = "../data_input/exclusions_by_missing_rate_proteins_individuals.RData")