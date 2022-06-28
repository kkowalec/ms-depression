###
## PRS analysis of MS MDD overlap
###

library(broom)
library(tidyverse)
library(fs)
library(glue)
z_path <- "/cifs/Z/arvhar/Project/MS_MDD/"



prs <- read_tsv("raw_data/transformed/joined_pgs.tsv")
ms <- read_tsv("raw_data/transformed/ms_cases.id") %>%
  mutate(ms = 1) %>% select(f.eid=1,ms)
ms <- left_join(ms, prs)

mdd_cidi_sf <- read_tsv("raw_data/transformed/mdd_cidi_sf.tsv")
mdd_icd <- read_tsv("raw_data/transformed/mdd_icd")
mdd_icd_or_cidi <- read_tsv("raw_data/transformed/mdd_icd_or_cidi.tsv")
mdd_sr <- read_tsv("raw_data/transformed/mdd_selfreported.tsv")
mdd_phq5 <- read_tsv("raw_data/transformed/mdd_phq_five.tsv")
mdd_phq10 <- read_tsv("raw_data/transformed/mdd_phq_ten.tsv")
mdd_phq_score <- read_tsv("raw_data/transformed/mdd_phq_score.tsv")
controls <- read_tsv("raw_data/transformed/controls.prs")
pass_qc <- read_tsv("raw_data/transformed/qc_cohort.tsv")
autoimmune <- read_tsv("raw_data/transformed/autoimmune_cases.tsv")



prs_analysis <- function(cases, ms, logistic = TRUE){
  
  if(logistic){
    
  data <- ms %>% 
    mutate(case = if_else(f.eid %in% cases$f.eid, 1,0))
  
  case_control <- data %>% 
    count(case)
  
  regression <- 
    data %>% 
    glm(case ~ mdd_no23noUKB_2018 +bmi+ sex + age + 
          PC1 +PC2 + PC3 + PC4 + PC5, data = ., family = "binomial") %>% 
    tidy(conf.int = TRUE, exponentiate = TRUE)
 
   regression %>% 
    mutate(ncase = case_control[[2,2]],
         ncontrol = case_control[[1,2]])
  
  } else {
    
    data <- 
      ms %>% 
      left_join(cases) %>% 
      filter(!is.na(tot_score))
    
    case <- data %>% 
      nrow()
    
    regression <- data %>% 
      lm(tot_score ~ mdd_no23noUKB_2018 + bmi + sex + age + 
            PC1 +PC2 + PC3 + PC4 + PC5, data = .) %>% 
      tidy(conf.int = TRUE)
    
    regression %>% 
      mutate(ncase = case,
             ncontrol = 0)
    
  }
  
}

aim1a_sex_stratified <- function(cases,ms, logistic = TRUE){
  if(logistic){
    
    ## Male
    data <- ms %>% 
      mutate(case = if_else(f.eid %in% cases$f.eid, 1,0)) %>% 
      filter(sex == 1)
    
    case_control <- data %>% 
      count(case)
    
    regression <- 
      data %>% 
      glm(case ~ mdd_no23noUKB_2018 +bmi + age + 
            PC1 +PC2 + PC3 + PC4 + PC5, data = ., family = "binomial") %>% 
      tidy(conf.int = TRUE, exponentiate = TRUE)
    
    male <- regression %>% 
      mutate(ncase = case_control[[2,2]],
             ncontrol = case_control[[1,2]]) %>% 
      mutate(sex = 'male')
    
    ## female
    #################################################################
    data <- ms %>% 
      mutate(case = if_else(f.eid %in% cases$f.eid, 1,0)) %>% 
      filter(sex == 2)
    
    case_control <- data %>% 
      count(case)
    
    regression <- 
      data %>% 
      glm(case ~ mdd_no23noUKB_2018 +bmi + age + 
            PC1 +PC2 + PC3 + PC4 + PC5, data = ., family = "binomial") %>% 
      tidy(conf.int = TRUE, exponentiate = TRUE)
    
    female <- regression %>% 
      mutate(ncase = case_control[[2,2]],
             ncontrol = case_control[[1,2]]) %>% 
      mutate(sex = 'female')
    
    bind_rows(male, female)
    
  } else {
  
  # Male
  data <- 
    ms %>% 
    filter(sex == 1)
    left_join(cases) %>% 
    filter(!is.na(tot_score))
  
  case <- data %>% 
    left_join(cases) %>% 
    filter(!is.na(tot_score)) %>% 
    nrow()
  
  regression <- data %>% 
    lm(tot_score ~ mdd_no23noUKB_2018 + bmi + sex + age + 
         PC1 +PC2 + PC3 + PC4 + PC5, data = .) %>% 
    tidy(conf.int = TRUE)
  
  male <- regression %>% 
    mutate(ncase = case,
           ncontrol = 0) %>% 
    mutate(sex = 'male')
  
  ## Female
  ################################################
  data <- 
    ms %>% 
    filter(sex == 2)
  left_join(cases) %>% 
    filter(!is.na(tot_score))
  
  case <- data %>% 
    left_join(cases) %>% 
    filter(!is.na(tot_score)) %>% 
    nrow()
  
  regression <- data %>% 
    lm(tot_score ~ mdd_no23noUKB_2018 + bmi + sex + age + 
         PC1 +PC2 + PC3 + PC4 + PC5, data = .) %>% 
    tidy(conf.int = TRUE)
  
  female <- regression %>% 
    mutate(ncase = case,
           ncontrol = 0) %>% 
    mutate(sex = 'female')
  
  bind_rows(male, female)
  
  
  }
}

cidi_sf <- prs_analysis(mdd_cidi_sf, ms)
icd <- prs_analysis(mdd_icd, ms)
icd_or_cidi <- prs_analysis(mdd_icd_or_cidi, ms)
selfreported <- prs_analysis(mdd_sr, ms)


ms2 <- ms %>% semi_join(mdd_phq_score)
phq5 <- prs_analysis(mdd_phq5,ms2)
phq10 <- prs_analysis(mdd_phq10,ms2)
phq_score <- prs_analysis(mdd_phq_score,ms, logistic = FALSE)


aim1_results <- list(cidi_sf, icd, icd_or_cidi, selfreported, phq5, phq10, phq_score)
definition <- c("cidi_sf", "icd", "icd_or_cidi", "selfreported", "phq5", "phq10", "phq_score")
nf <- function(tbl,name){tbl %>% 
    mutate(name = {{name}}) %>% 
    filter(term %in% c('mdd_no23noUKB_2018', "bmi")) %>% 
    select(estimate, conf.low, conf.high,ncase, ncontrol, p.value,term, name)
  }

aim1_clean_results <- map2(aim1_results, definition,nf) %>% 
  reduce(bind_rows)


## sex stratified models
cidi_sf <- aim1a_sex_stratified(mdd_cidi_sf, ms)


phq_stats <- 
  mdd_phq_score %>% 
  semi_join(ms) %>% 
  summarise(mean = mean(tot_score), sd = sd(tot_score))

mdd_results <- bind_rows(c,a, b,d,f,g,h)


# Aim 1B

aim_1b <- function(ms, mdd, autoimmune){
  
  
  case <- mdd %>% 
    mutate(case = if_else(f.eid %in% ms$f.eid,1,0)) %>%
    filter(case == 1) %>%
    select(1, case)
  
  
  control <-  mdd %>%
    anti_join(case) %>% 
    anti_join(autoimmune) %>% 
    mutate(case = 0) %>%  
    select(1, case)
  
  dataset <- bind_rows(case, control)
  
  
  case_control <- dataset %>%
    count(case)
  
  regression <- dataset %>%
    left_join(prs) %>% 
    glm(case ~ mdd_no23noUKB_2018 + bmi + PC1 + PC2 + PC3 + PC4 + PC5 + age + sex, family = "binomial", data = .) %>%
    tidy(exponentiate = TRUE, conf.int = TRUE)
  
  regression %>%
    mutate(ncase = case_control[[2,2]],
           ncontrol = case_control[[1,2]]
    )
}

aim_1c <- function(ms, controls, mdd){
 controls <- controls %>% mutate(case = 0) %>% select(f.eid, case)
 
 ms_mdd <- ms %>% 
   semi_join(mdd) %>% 
   mutate(case = 1) %>% 
   select(f.eid, case)
 
 dataset <- bind_rows(controls, ms_mdd)
 
 case_control <- dataset %>% count(case)
 
 regression <- dataset %>%
   left_join(prs) %>% 
   glm(case ~ mdd_no23noUKB_2018 + bmi + PC1 + PC2 + PC3 + PC4 + PC5 + age + sex, family = "binomial", data = .) %>% 
   tidy(exponentiate = TRUE, conf.int = TRUE)
 
 
 regression %>%
   mutate(ncase = case_control[[2,2]],
          ncontrol = case_control[[1,2]]
   )
 
 

}

md_definitions <- list(mdd_cidi_sf, mdd_icd, mdd_icd_or_cidi, mdd_sr ,mdd_phq5,
                       mdd_phq10, mdd_phq_score )


aim_1b_results <- pmap(
  list(
    list(ms),
    md_definitions,
    list(autoimmune)
  ),
  aim_1b)


aim_1c_results <- pmap(
  list(
    list(ms),
    list(controls), 
    md_definitions
    ),
  aim_1c)


aim_1a_sexstrat <- pmap(
  list(
    md_definitions,
    list(ms)
  ),
  aim1a_sex_stratified)

aim_1_results <- list(aim_1b_results, aim_1c_results, aim_1a_sexstrat)
save(aim_1_results, file = glue(z_path,"results/aim1.RData"))


#
mdd_cidi_sf %>% 
  anti_join(autoimmune) %>% 
  left_join(mdd_phq_score) %>% 
  filter(!is.na(tot_score)) %>% 
  left_join(prs) %>% 
  lm(tot_score ~ mdd_no23noUKB_2018  + bmi +  PC1 + PC2 + PC3 + PC4 + PC5 + age + sex, data = .) %>% 
  tidy() %>% 
  filter(term %in% c('mdd_no23noUKB_2018', 'bmi')) %>% 
  write_tsv(glue(z_path, 'results/mdd_phq9{Sys.Date()}.tsv'))


  

