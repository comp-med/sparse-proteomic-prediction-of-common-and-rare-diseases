###################################################
#### Script Name: generate_input_file
#### Author: Julia Carrasco Zanini Sanchez
#### Description: generates inpute file to run ML framework in a job array 
#### Date : 20/07/2022

rm(list = ls())

library(tidyr)
library(dplyr)

## set working directory 
setwd("/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/01_inc_dz_phenos/bin/") 

## read in disease list 
dz.list <- read.delim("../data_input/atlas308_subc_inc6months_inc_dz_counts_3_10yrs.txt", sep = "\t")
head(dz.list)

dz.input <- dz.list %>% 
  select(inc_dz_name, dz_date) %>% 
  mutate(inc_yrs=rep(10,nrow(dz.list)))

# dz.input <- rbind(dz.input, dz.list %>% 
#   select(inc_dz_name, dz_date) %>% 
#   mutate(inc_yrs=rep(5,nrow(dz.list)))) %>% 
#   filter(!inc_dz_name%in%c("atlas_inc_ms","atlas_inc_sle"))

write.table(dz.input, file = "../../02_inc_dz_prediction/data_input/atlas308_file_input_params.txt", sep="\t", col.names = F, row.names = F, quote = F)

## diseases with fewer cases include consortium selected cases
dz.list.cc <- read.delim("../data_input/atlas308_consort_inc6months_inc_dz_counts_3_10yrs.txt", sep = "\t")

dz.input.cc <- dz.list.cc %>% 
  select(inc_dz_name, dz_date) %>% 
  mutate(inc_yrs=rep(10,nrow(dz.list.cc)))

dz.input.cc <- rbind(dz.input.cc, dz.list.cc %>% 
        select(inc_dz_name, dz_date) %>% 
        mutate(inc_yrs=rep(5,nrow(dz.list.cc))))
write.table(dz.input.cc, file = "../../02_inc_dz_prediction/data_input/atlas308_file_input_params_consort.txt", sep="\t", col.names = F, row.names = F, quote = F)

# list for PRS prediction models
prs.input <- dz.input %>%
  filter(inc_dz_name%in%c("atlas_inc_asthma","atlas_inc_af",
                          "atlas_inc_pri_bowel","atlas_inc_pri_breast",
                          "atlas_inc_chd_nos","atlas_inc_crohns",
                          "atlas_inc_pri_ovarian", "atlas_inc_fracture_hip",
                          "atlas_inc_fracture_wrist","atlas_inc_isch_stroke",
                          "atlas_inc_pri_skin",
                          "atlas_inc_oa","atlas_inc_parkinsons",
                          "atlas_inc_glaucoma","atlas_inc_pri_prost",
                          "atlas_inc_psoriasis","atlas_inc_rha",
                          "atlas_inc_diabetes_t2","atlas_inc_ulc_colitis",
                          "atlas_inc_vte_ex_pe"))
write.table(prs.input, file = "../../02_inc_dz_prediction/data_input/prs_input_params.txt", sep="\t", col.names = F, row.names = F, quote = F)
