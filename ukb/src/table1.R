library(broom)
library(tidyverse)
library(fs)
library(glue)
z_path <- "/cifs/Z/arvhar/Project/MS_MDD/"

prs <- read_tsv("raw_data/transformed/joined_pgs.tsv")
ms <- read_tsv("raw_data/transformed/ms_cases.id") %>% select(f.eid=IID)
mdd_icd_or_cidi <- read_tsv("raw_data/transformed/mdd_icd_or_cidi.tsv")
controls <- read_tsv("raw_data/transformed/controls.prs")
ms <- read_tsv("raw_data/transformed/ms_cases.id") %>% rename(f.eid = IID)
mdd_phq_score <- read_tsv("raw_data/transformed/mdd_phq_score.tsv")
mdd_sr <- read_tsv("raw_data/transformed/mdd_selfreported.tsv")
mdd_phq5 <- read_tsv("raw_data/transformed/mdd_phq_five.tsv")
mdd_phq10 <- read_tsv("raw_data/transformed/mdd_phq_ten.tsv")

# lifesyle variables
inc_education <- read_tsv("raw_data/ukb_edu_years.tsv") %>% tibble() %>% 
  select(-household_income)
income <- read_tsv("raw_data/household_income.tsv")
bmi <-  read_tsv("raw_data/bmi.tsv") %>% select(f.eid, bmi_measure=bmi)
eversmoke <- read_tsv('raw_data/transformed/eversmoked.tsv')

# males = 1
count(ms, sex)

# Coverage with PHQ-9
mdd_icd_or_cidi %>% 
  left_join(mdd_phq_score) %>% 
  summarise(sum(!is.na(tot_score)))

ms %>% 
  left_join(mdd_phq_score) %>% 
  summarise(sum(!is.na(tot_score)))

controls %>% 
  left_join(mdd_phq_score) %>% 
  summarise(sum(!is.na(tot_score)))

## Table S1
data <- controls
name <- 'controls'
inc_edu <- function(data,name){
  income <- data %>% 
    left_join(income) %>% 
    count(household_income) %>% 
    rename({{name}} := n, measure = household_income) %>% 
    mutate(measure = if_else(is.na(measure), 'missing_inc', measure))
  
  edu <- data %>% 
    left_join(inc_education) %>% 
    count(education_level) %>% 
    rename({{name}} := n, measure = education_level) %>% 
    mutate(measure = if_else(is.na(measure), 'missing_edu', measure))
  bind_rows(income,edu)
} 





make_table <- function(data,name){
  data %>% 
    left_join(mdd_phq_score) %>% 
    left_join(prs) %>% 
    left_join(bmi) %>% 
    left_join(inc_education) %>% 
    summarise(female = sum(sex==2, na.rm=T),
            pct_female = female/nrow(.),
            mean_age  = mean(age,na.rm=T),
            sd_ag =  sd(age,na.rm=T),
            lifetime_md = sum(f.eid %in% mdd_icd_or_cidi$f.eid),
            selfreported_md = sum(f.eid %in% mdd_sr$f.eid),
            ever_smoker = sum(f.eid %in% eversmoke$f.eid),
            pct_ever_smoker = ever_smoker/nrow(.),
            mean_eyears= mean(edu_years,na.rm=T),sd_eyears= sd(edu_years,na.rm=T),
            mean_bmi = mean(bmi_measure,na.rm =T), sd_bmi = sd(bmi_measure,na.rm=T),
            mean_phq = mean(tot_score, na.rm=T),
            sd_phq = sd(tot_score, na.rm=T),
            md_phq5 = sum(f.eid %in% mdd_phq5$f.eid),
            md_phq10 = sum(f.eid %in% mdd_phq10$f.eid),
            md_prs_mean = mean(mdd_no23noUKB_2018,na.rm=T), md_prs_sd = sd(mdd_no23noUKB_2018,na.rm=T),
            bmi_prs_mean = mean(bmi,na.rm=T), bmi_prs_sd = sd(bmi, na.rm=T),
            
  ) %>% 
    pivot_longer(everything()) %>% rename({{ name }} := value)
}


phenos <- list(ms,mdd_icd_or_cidi, controls)
names <- c('ms','mdd','controls')

table1 <- map2(phenos, names, make_table) %>% 
  reduce(inner_join, by = "name")

table1_inc_edu <- map2(phenos, names, inc_edu) %>% 
  reduce(inner_join, by = "measure")

write_tsv(table1_inc_edu,'results/table1_inc_edu.tsv')
write_tsv(table1,'results/table1.tsv')
