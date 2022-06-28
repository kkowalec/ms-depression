library(tidyverse)
library(data.table)
library(fs)
source("functions/functions.R")
load("produced_data/exclusion_criteria_for_controls.RData")
pass_qc <- read_tsv("raw_data/transformed/qc_cohort.tsv")

# individuals meeting criteria for anxiety or MDD

all_mdd <- read_tsv("raw_data/transformed/all_mdd.tsv")
all_anx <- read_tsv("raw_data/transformed/all_anx.tsv")
mdd_or_anx <- full_join(all_mdd, all_anx)

# Get phenotype data

icd_data <- extract_ukb('41270')
meds_data <- extract_ukb('20003')
selfrep_data <- extract_ukb('20002')
selfrep_cancer <- extract_ukb('2453')

#################
## Derive all individuals with ICD codes to exclude from controls
##################

icd_phenotypes <- 
  map2(joined_definitions$icd, 
       joined_definitions$phenotype, .f = check_for_code, data = icd_data) %>% 
  reduce(full_join)

names(icd_phenotypes) <- paste0(colnames(icd_phenotypes), "_icd")

icd_phenotypes <- 
  icd_phenotypes %>% 
  rename(f.eid = f.eid_icd)

# Derive all cases with self-reported code

selfrep_phenotypes <- 
  map2(joined_definitions$selfrep, 
       joined_definitions$phenotype, .f = check_for_code, data = selfrep_data) %>% 
  reduce(full_join)

names(selfrep_phenotypes) <- paste0(colnames(selfrep_phenotypes), "_sr")

selfrep_phenotypes <- 
  selfrep_phenotypes %>% 
  rename(f.eid = f.eid_sr)


# Join dataframes together

full_exclusion <- full_join(selfrep_phenotypes, icd_phenotypes, by ="f.eid")

# Define ever-diagnosed with cancer on any occasion (some individuals have multiple visits)

selfrep_cancer <- 
  selfrep_cancer %>%
  mutate(across(2:4, ~replace_na(.,0))) %>% 
  filter(if_any(-1, ~ . == 1)) %>% 
  mutate(cancer_sr_temp = 1) %>% 
  select(1, cancer_sr_temp)

# Update cancer_sr column and cancer_both

exclusion <- 
  full_exclusion %>% 
  full_join(., selfrep_cancer) %>% 
  mutate(across(-1, ~replace_na(., 0))) %>% 
  mutate(cancer_sr = case_when(cancer_sr_temp == 1 ~ 1,
                               cancer_sr_temp != 1 ~ 0))

# Define controls as individuals in cohort WITHOUT any of the ICD codes and
# selfreported disorders extracted in the start of the scri

controls <- 
  pass_qc %>% 
  anti_join(.,mdd_or_anx, by = "f.eid") %>% 
  anti_join(., exclusion, by = "f.eid") %>% 
  filter(in.kinship.table == 0) %>% 
  mutate(controls = 1) %>% 
  select(f.eid, controls)

write_tsv(controls, "raw_data/transformed/controls.prs")