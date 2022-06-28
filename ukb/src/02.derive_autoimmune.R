library(tidyverse)
library(data.table)
source("functions/functions.R")

qc_cohort <- read_tsv("raw_data/transformed/qc_cohort.tsv")
load("produced_data/exclusion_criteria_for_controls.RData")
load("produced_data/immune_disorder_meds.RData")

# see functions/functions.R for these helper functions.
icd_data <- extract_ukb('41270')
meds_data <- extract_ukb('20003')
selfrep_data <- extract_ukb('20002')
relatedness <- fread("raw_data/ukb_related.dat") %>% tibble()


ms <- 
  derive_ukb_n(icd_codes = "G35", selfrep_codes = "1261", 
               medicine_codes = as.integer(auto_immune_meds$ms), col_name = "multiple_sclerosis") 

ibd <- 
  derive_ukb_n(icd_codes = c("K50", "K51"), selfrep_codes = c("1463","1462"), 
               medicine_codes = auto_immune_meds$ibd, col_name ="inflammatory_bowel_disorder") 

ra <- 
  derive_ukb_n(icd_codes = c("M05", "M06"), selfrep_codes = "1464", 
               medicine_codes = auto_immune_meds$ra, col_name ="rheumatoid_arthritis") 




# The probable_case function is used to only retain cases that meet atleast 2/3 criteria.

ra %>% 
  probable_case() %>%
  mutate(fid = f.eid) %>% 
  select(f.eid, fid) %>% 
  write_tsv(., "raw_data/transformed/ra_cases.id", col_names = FALSE)

ms %>% 
  probable_case() %>% 
  mutate(fid = f.eid) %>% 
  select(f.eid, fid) %>% 
  write_tsv(., "raw_data/transformed/ms_cases.id", col_names = FALSE)

ibd %>% 
  probable_case() %>% 
  mutate(fid = f.eid) %>% 
  select(f.eid, fid) %>% 
  write_tsv(., "raw_data/transformed/ibd_cases.id", col_names = FALSE)

# identify any individuals with an autoimmune disorderr

icd_codes <- filter(joined_definitions, phenotype == "autoimmune_disorder")
sr <- filter(joined_definitions, phenotype == "autoimmune_disorder")

emr <- check_for_code(codes = icd_codes$icd[[1]], col_name = "autoimmune", icd_data)
selfrep <- check_for_code(codes = sr$selfrep[[1]], col_name = "autoimmune", selfrep_data)

any_autoimmune <- full_join(emr, selfrep, by = "f.eid")

# save list of ids

write_tsv(any_autoimmune, "raw_data/transformed/autoimmune_cases.tsv")

# for MS, remove one individual in each related pair

ms <-  read_tsv("raw_data/transformed/ms_cases.id", col_names = FALSE)

# find pairs of related individuals where one has MS
relatedness %>% 
  filter(ID1 %in% ms$X1 & ID2 %in% ms$X1)
