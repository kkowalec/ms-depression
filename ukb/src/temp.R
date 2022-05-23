library(tidyverse)
library(data.table)
source("functions/functions.R")
library(broom)

# AIM:
# correlation between BMI pgs and BMI in four subgroups
#
mdd_icd_or_cidi <- read_tsv("raw_data/transformed/mdd_icd_or_cidi.tsv")
prs <- read_tsv("raw_data/transformed/joined_pgs.tsv")
ms <- read_tsv("raw_data/transformed/ms_cases.id") %>%
  mutate(ms = 1) %>% select(f.eid=1,ms)

bmi <- extract_ukb("21001") %>% select(f.eid,f.21001.0.0)

ms_dep <- semi_join(ms,mdd_icd_or_cidi)
dep_no_ms <- anti_join(mdd_icd_or_cidi,ms)




bmi_cor <- function(tbl) {
  df <- tbl %>% 
    left_join(bmi) %>% 
    left_join(prs) %>% 
    filter(!is.na(bmi) & !is.na(f.21001.0.0))
  
  tibble(cor = stats::cor(df$bmi, df$f.21001.0.0))
  
}

names <- c("mdd", "healthy", "ms_dep", "dep_no_ms")
correlations <- list(mdd_icd_or_cidi,controls,ms_dep,dep_no_ms) %>% 
  map_df(bmi_cor) %>% 
  add_column(pheno = names)

colnames(prs)

# What is the r2 for bmi and MDD prs in mdd and ms vs mdd no ms?

pheno <- mutate(ms,case = if_else(f.eid %in% mdd_icd_or_cidi$f.eid,1,0)) %>% 
  left_join(prs)

# model with mdd prs
ll <- glm(case ~ mdd_no23noUKB_2018 + PC1 + PC2 + PC3 + PC4 + PC5 + sex + age, family = "binomial", data = pheno) %>% 
  glance() %>% pull(logLik)

# model with bmi prs
ll2 <- glm(case ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex + age, family = "binomial", data = pheno) %>% 
  glance() %>% pull(logLik)

mdd_r2 <- 1-(ll/ll2)


ll <- glm(case ~ bmi + PC1 + PC2 + PC3 + PC4 + PC5 + sex + age, family = "binomial", data = pheno) %>% 
  glance() %>% pull(logLik)

# model with bmi prs
ll2 <- glm(case ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex + age, family = "binomial", data = pheno) %>% 
  glance() %>% pull(logLik)

bmi_r2 <- 1-(ll/ll2)

