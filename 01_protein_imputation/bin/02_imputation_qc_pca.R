###################################################
#### Script Name: 02_imputation_qc_pca
#### Author: Julia Carrasco Zanini Sanchez
#### Description: script to do a PCA on imputed npx protein values 
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

## set up working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/00_protein_imputation/bin/") 

## load imputed NPX values
load("../data_input/Imputed_NPX_missForest_Cardiometabolic_II.RData")
imp.npx <- m.forest$ximp[,-c(1,2)]
for (i in c("Inflammation_II","Neurology_II","Oncology_II","Cardiometabolic","Inflammation","Neurology","Oncology")) {
  load(paste0("../data_input/Imputed_NPX_missForest_",i,".RData"))
  imp.npx <- cbind(imp.npx, m.forest$ximp[,-c(1,2)])
}

## merge to age and sex
ids <- read.delim("../data_input/temp_imputation_ids.txt", sep = "\t")
ol.npx <- cbind(ids, imp.npx)

## run pca 
pca.prots <- prcomp(ol.npx[,5:ncol(ol.npx)])
library(ggfortify)

## load exclusions by missing rate
load("../data_input/exclusions_by_missing_rate_proteins_individuals.RData")

ol.npx$exclusion.50 <- as.factor(ifelse(ol.npx$eid_20361%in%ex.mis.rate$individuals[which(ex.mis.rate$individuals.mis.rate>50)],1,0))

## check whether imputation introduced any segregation by missingness -- no evidence of this
autoplot(pca.prots, data = ol.npx, colour = "mis.60")

write.table(ol.npx, file = "../data_input/Imputed_1536_Expansion_NPX_proteins.txt", sep = "\t", row.names = F)
