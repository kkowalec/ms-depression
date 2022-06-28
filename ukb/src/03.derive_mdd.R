library(data.table)
library(fs)
library(glue)
library(tidyverse)

source("functions/functions.R")
source("functions/repro_functions.R")


###
## READ IN DATA

cidi_cols <- c("20446", "20441", "20436", "20439", "20440", "20449", "20536", "20532", "20435", "20450", "20437")
phq9_cols <- c(20514, 20510, 20517,  20519, 20511, 20507, 20508, 20518, 20513)

sf_data <- extract_ukb(cidi_cols)
phq_data <- extract_ukb(as.character(phq9_cols))

###
## APPLY DEFINITIONS FROM DAVIS et al, 2020

cidi_sf_cases <- 
  sf_data %>% 
  # rename columns to make it easier to work with
  rename_with(~str_sub(.x, start = 1, end = 7)) %>% 
  filter(f.20446 == 1 | f.20441 == 1) %>% 
  filter(f.20436 > 2 & f.20439 > 2 & f.20440 > 2) %>% 
  select(f.eid, f.20446, f.20441, f.20449, f.20536, f.20532, f.20435, f.20450, f.20437) %>% 
  # Prefer not to answer and one other answer option is coded as -121 or -818. Remove these.
  mutate(across(-f.eid, ~ifelse(.x == -121 | .x == -818, 0, .x))) %>% 
  mutate(tot_score = rowSums(.[-1])) %>% 
  filter(tot_score >= 5) %>% 
  mutate(md_cidi= 1) %>% 
  select(f.eid, md_cidi)


phq_totscore <- 
  phq_data %>% 
  rename_with(~str_sub(.x, start = 1, end = 7)) %>%
  # Remove those with NA
  filter(if_all(everything(), ~!is.na(.x))) %>% 
  # remove prefer not to answer, or missing
  mutate(across(-f.eid, ~ifelse(.x == -121 | .x == -818, 1, .x))) %>% 
  # PHQ-9 seems coded differetly from others, when using cut-off 5 and 10.
  # normalise to go for
  mutate(across(-f.eid, ~.x-1)) %>% 
  mutate(tot_score = rowSums(.[-1], na.rm = TRUE)) %>% 
  select(f.eid, tot_score)

# individuals with tot_score >= 5

phq_five <-  
  phq_totscore %>% filter(tot_score >= 5)

# individuals with tot_score >= 10

phq_10 <- 
  phq_totscore %>% filter(tot_score >= 10)


#######
## WRITE OUT RESULTS

write_tsv(cidi_sf_cases,"raw_data/transformed/mdd_cidi_sf.tsv")
write_tsv(phq_five,"raw_data/transformed/mdd_phq_five.tsv")
write_tsv(phq_10,"raw_data/transformed/mdd_phq_ten.tsv")
write_tsv(phq_totscore,"raw_data/transformed/mdd_phq_score.tsv")


#######
## Define MDD based on ICD codes

icd_data <- extract_ukb('41270')
icd_phenotypes <- check_for_code(codes = c("F32", "F33"), col_name = "icd_mdd", data =  icd_data)

# save the results

write_tsv(icd_phenotypes,"raw_data/transformed/mdd_icd")


## ICD or CIDI-sf (merged)

full_join(icd_phenotypes, cidi_sf_cases) %>% 
  mutate(icd_or_cidi = 1) %>% 
  select(f.eid, icd_or_cidi) %>% 
  write_tsv("raw_data/transformed/mdd_icd_or_cidi.tsv")


#####
## Selfreported MDD (MHQ and baseline)

mhq_sr <- extract_ukb('20544')
md_mhq_selfrep <- check_for_code(codes = "11", data = mhq_sr, col_name = "mhq_selfrep") 

selfrep_data <- extract_ukb('20002')
selfrep_baseline <- check_for_code(codes = "1286", data = selfrep_data, col_name = "baseline_selfrep")



# Cardinal symptoms, from CIDI-sf

cardinal_symptoms <- 
  sf_data %>% 
  filter(if_any(c(f.20441.0.0,f.20446.0.0), ~.== 1)) %>% 
  mutate(cardinal = 1) %>% 
  select(f.eid, cardinal)

# Read in datafield 20216 - Smith et al definition of MDD

mdd_touchscreen <- extract_ukb('20126') %>% 
  filter(f.20126.0.0 %in% c(3,4,5)) %>% 
  rename(touchscreen =f.20126.0.0)


all_mdd_definitions <- list(selfrep_baseline, md_mhq_selfrep,
                            icd_phenotypes, cidi_sf_cases, phq_totscore,
                            cardinal_symptoms,mdd_touchscreen)

# Bind together into one dataframe

df_mdd <- reduce(all_mdd_definitions, full_join, by = "f.eid")

joined_final <- df_mdd %>% select(-tot_score) %>% filter(if_any(-1, ~!is.na(.)))

# save definitions

write_tsv(joined_final,"raw_data/transformed/all_mdd.tsv")

# How many MDD cases from each definition?

df_mdd %>% 
  select(-f.eid,-tot_score) %>% 
  filter(if_any(everything(), ~!is.na(.))) %>% 
  map_df(~sum(!is.na(.))) %>% 
  write_tsv("produced_data/mdd_case_n.tsv")
  




  