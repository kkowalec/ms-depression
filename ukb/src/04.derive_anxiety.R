library(tidyverse)
library(data.table)
library(fs)
library(glue)

source("functions/functions.R")
source("functions/repro_functions.R")

qc_cohort <- read_tsv("raw_data/transformed/qc_cohort.tsv")
load("produced_data/immune_disorder_meds.RData")


###
## CIDI-SF Columns

anx_cols <- c("20421", "20420", "20538", "20425", "20542","20543", "20540", "20541", "20539", "20537", "20418",
              "20426", "20423", "20429", "20419", "20422", "20417", "20427") 

cidi_sf_anx <- extract_ukb(anx_cols)

# DAVIS 2020 definition

anx_cidi <- 
  cidi_sf_anx %>% 
  semi_join(qc_cohort) %>% 
  filter(f.20421.0.0 == 1) %>% 
  filter(f.20420.0.0 >= 6 | f.20420.0.0 == -999 ) %>% 
  filter(f.20538.0.0 == 1) %>% 
  #filter(f.20425.0.0 == 1 | f.20542.0.0 == 1) %>% 
  filter(f.20543.0.0 == 2 | f.20540.0.0 == 1) %>% 
  filter(f.20541.0.0 == 1 | f.20539.0.0 == 3 | f.20537.0.0 == 3) %>% 
  filter(f.20418.0.0 > 1) %>% 
  mutate(across(c("f.20426.0.0", "f.20423.0.0", "f.20429.0.0", 
                  "f.20419.0.0", "f.20422.0.0", "f.20417.0.0", "f.20427.0.0"), ~if_else(.x == 1, 1, 0))) %>% 
  rowwise() %>% 
  mutate(n_symptoms = sum(
    c_across(c("f.20426.0.0", "f.20423.0.0", "f.20429.0.0", 
               "f.20419.0.0", "f.20422.0.0", "f.20417.0.0", "f.20427.0.0")))) %>% 
  ungroup() %>% 
  filter(n_symptoms >= 3) %>% 
  mutate(anx_cidi = 1) %>% 
  select(f.eid, anx_cidi)

# save individuals

write_tsv(anx_cidi,"raw_data/transformed/anx_cidi_sf.tsv")


###
## GAD-7

gad7_cols <-c(20506, 20509, 20520, 20515, 20516, 20505, 20512) 

gad7 <- extract_ukb(as.character(gad7_cols))

gad7_score <- 
  gad7 %>% 
  # Remove those with NA
  filter(if_all(everything(), ~!is.na(.x))) %>% 
  # remove prefer not to answer, or missing
  mutate(across(-f.eid, ~ifelse(.x == -121 | .x == -818, 1, .x))) %>% 
  mutate(across(-f.eid, ~.x-1)) %>% 
  mutate(tot_score = rowSums(.[-1], na.rm = TRUE),
         tot_score = tot_score) %>% 
  select(f.eid, tot_score)


# save results

write_tsv(gad7_score,"raw_data/transformed/anx_gad7_score.tsv")

gad7_score %>% 
  filter(tot_score >= 7) %>% 
  write_tsv("raw_data/transformed/anx_gad7_7.tsv")

gad7_score %>% 
  filter(tot_score >= 10) %>% 
  write_tsv("raw_data/transformed/anx_gad7_10.tsv")

gad7 %>% 
  mutate(across(-f.eid, ~ifelse(.x == -121 | .x == -818, 1, .x))) %>% 
  mutate(across(-f.eid, ~.x-1)) %>% 
  mutate(tot_score = f.20506.0.0 + f.20509.0.0) %>% 
  filter(!is.na(tot_score)) %>% 
  select(f.eid, tot_score) %>% 
  write_tsv("raw_data/transformed/anx_gad7_firstwo.tsv")


######
## ICD based definitions

icd_data <- extract_ukb('41270')

check_for_code(c("F411", "F40", "F410"), "dsm_v", icd_data) %>% 
  full_join(anx_cidi, by = "f.eid") %>% 
  write_tsv("raw_data/transformed/anx_def2.tsv")

check_for_code(c("F411", "F40", "F410", "F42", "F431"), "dsm_iv", icd_data) %>% 
  full_join(anx_cidi, by = "f.eid") %>% 
  write_tsv("raw_data/transformed/anx_def3.tsv")

icd <- check_for_code(c("F411", "F40", "F410", "F42", "F431", "F419"), "icd", icd_data)

full_join(icd,anx_cidi, by = "f.eid") %>% 
  write_tsv("raw_data/transformed/anx_def4.tsv")


######
## self-reported definitions from MHQ and baseline

mhq_sr <- extract_ukb('20544')
baseline_sr <- extract_ukb('20002')


# define self-reported anxiety using code 15, merge with baseline selfreport

selfrep <- mhq_sr %>% 
  filter(if_any(-1, ~. == 15)) %>% 
  mutate(mhq_sr_anx = 1) %>% 
  select(1, mhq_sr_anx) %>% 
  full_join(check_for_code( "1287","anx",data = baseline_sr))

# save results

write_tsv(selfrep, "raw_data/transformed/anx_selfreported.tsv")
combined_anx <- full_join(selfrep,icd) %>% full_join(anx_cidi)
combined_anx %>%  write_tsv("raw_data/transformed/all_anx.tsv")


list(selfrep, anx_cidi,gad7_score,icd) %>% 
  reduce(full_join, by = "f.eid") %>% 
  mutate(gad7 = ifelse(tot_score >= 7, 1,NA),
         gad10 = ifelse(tot_score >= 10, 1,NA)) %>% 
  select(-f.eid, -tot_score) %>% 
  map_df(~sum(!is.na(.))) %>% 
  write_tsv("produced_data/anxiety_case_n.tsv")
  