library(tidyverse)
library(data.table)
library(fs)


###
## Define a subset of full UKB that pass QC filters.
## European, pass genetic QC, has not withdrawn their consent
###


withdrawn <- fread("raw_data/withdrawal_list_20210809.csv") %>% tibble() %>% rename(iid = V1)


# Non-european participants
ancestry_outlier <- 
  read_delim( "raw_data/transformed/ancestraloutliers_3sd.fid.iid", delim = " ", col_names = FALSE) %>% 
  rename(fid = X1, iid = X2) 


metadata <- 
  fread("raw_data/transformed/ukb2222_cal_v2_s488364_w_header_w_sqc_v2.txt") %>% 
  tibble() %>% 
  rename(iid = IID)


qced_cohort <- 
  metadata %>% 
  # Remove individuals who have withdrawn their participation
  anti_join(., withdrawn, by = "iid") %>% 
  # Did not pass genetic QC
  filter(het.missing.outliers == 0) %>% 
  # Did not pass genetic QC
  filter(sample.qc.missing.rate < 0.02) %>% 
  # Genetically inferred gender should be same as reported gender
  filter(Submitted.Gender == Inferred.Gender) %>% 
  # keep only individuals of european ancestry
  filter(!(iid %in% ancestry_outlier$iid)) %>% 
  rename(f.eid = iid) %>% 
  filter(SEX != 0) %>% 
  select(f.eid, in.kinship.table, SEX)


write_tsv(qced_cohort, "raw_data/transformed/qc_cohort.tsv")