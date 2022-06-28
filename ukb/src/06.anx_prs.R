
## PRS analysis of MS ANX overlap
###
library(broom)
library(tidyverse)
library(fs)


#
# READ IN PRS and MS data
#
prs <- read_tsv("raw_data/transformed/joined_pgs.tsv")
anx_prs <- read_tsv('raw_data/transformed/levey_anx_sr_banded.sscore') %>% 
  mutate(anx_2020_sr = as.numeric(scale(SCORE1_AVG))) %>% 
  select(f.eid = IID, anx_2020_sr)

prs <- left_join(prs, anx_prs)
ms <- read_tsv("raw_data/transformed/ms_cases.id") %>%
  mutate(ms = 1) %>% select(f.eid=1,ms)
ms <- left_join(ms, prs)

# Read in anxiety definitions
anx_cidi <- read_tsv("raw_data/transformed/anx_cidi_sf.tsv")
anx_def2 <- read_tsv("raw_data/transformed/anx_def2.tsv")
anx_def3 <- read_tsv("raw_data/transformed/anx_def3.tsv")
anx_def4 <- read_tsv("raw_data/transformed/anx_def4.tsv")
anx_sr <- read_tsv("raw_data/transformed/anx_selfreported.tsv")
anx_gad7 <- read_tsv("raw_data/transformed/anx_gad7_7.tsv")
anx_gad10 <- read_tsv("raw_data/transformed/anx_gad7_10.tsv")
anx_gad_f2 <- read_tsv("raw_data/transformed/anx_gad7_firstwo.tsv")
anx_gad_tot <- read_tsv("raw_data/transformed/anx_gad7_score.tsv")


# define regression formulas. Use two anxiety definitions for anx PRS
gad_score <- as.formula(case ~ anx_2020 + sex + age + PC1 +PC2 + PC3 + PC4 + PC5)
anx_selfrep <- as.formula(case ~ anx_2020_sr + sex + age + PC1 +PC2 + PC3 + PC4 + PC5)

# list the different definitions, and their names
basic_definitions <- list(anx_cidi,anx_def2, anx_def3, anx_def4, anx_sr)
basic_names <- c('cidi','def2', 'def3', 'def4', 'self_reported')

prs_analysis <- function(cases, ms, formula,logistic = TRUE){
  
  if(logistic){
    
    data <- ms %>% 
      mutate(case = if_else(f.eid %in% cases$f.eid, 1,0))
    
    k <- data %>% 
      count(case)
    
    l <- data %>% 
      glm(formula, data = ., family = "binomial") %>% 
      tidy(conf.int = TRUE, exponentiate = TRUE)
    list(k,l)
    
  } else {
    
    data <- 
      ms %>% 
      left_join(cases) %>% 
      filter(!is.na(tot_score))
    
    k <- tibble(case = 1, n = nrow(cases)) %>% 
      add_row(case = 0, n=0)
    
    l <- 
      data %>% 
      rename(case=tot_score) %>% 
      lm(formula, data = .) %>% 
      tidy(conf.int = TRUE)
    list(k,l)
    
  }
  
}
extract <- function(res, term, name){
  
  sample_size <- 
    res[[1]] %>% 
    pivot_wider(values_from = n, names_from = case) %>% 
    rename(control =1, case = 2)
  
  coefficient <- 
    res[[2]] %>% 
    filter(term == {{ term }} ) %>%  
    select(estimate, std.error, p.value, conf.low, conf.high)
  
  bind_cols(sample_size, coefficient) %>% add_column(definition = {{ name }}) %>% 
    select(definition, everything())
  
}

# run regressions over different definitions
res_gad <- map(basic_definitions, prs_analysis, ms=ms, formula=gad_score)
res_sr <- map(basic_definitions, prs_analysis, ms=ms, formula=anx_selfrep)

# extract relevant data
munged_gad <- map2_dfr(res_gad, basic_names, extract, term = 'anx_2020') %>% 
  add_column(prs = 'gad2')
munged_selrep <- map2_dfr(res_sr, basic_names, extract, term = 'anx_2020_sr') %>% 
  add_column(prs = 'self-report')

prs_res <- bind_rows(munged_gad, munged_selrep)
write_tsv(prs_res, glue(z_path, "results/anx_prs{Sys.Date()}.tsv"))






# Restrict MS to those with a GAD-7 score
ms2 <- ms %>% semi_join(anx_gad_f2, by = "f.eid")

gad7_7res <- prs_analysis(anx_gad7, ms, formula = gad_score)
gad7_10res <- prs_analysis(anx_gad10, ms,formula = gad_score)
gad7_firsttwo <- prs_analysis(anx_gad_f2, ms, logistic = FALSE, formula = gad_score)
gad7_totscore <- prs_analysis(anx_gad_tot, ms, logistic = FALSE, formula = gad_score)

gad7_names <- c('gad7', 'gad10','gad_f2', 'gad_tot')
gad7_res <- list(gad7_7res,gad7_10res, gad7_firsttwo,gad7_totscore)
munged_gad <- map2_dfr(gad7_res, gad7_names, extract, term = 'anx_2020')





###
## Average GAD7 and GAD7 score
###

stats1 <- 
  anx_gad_tot %>% 
  semi_join(ms) %>% 
  summarise(definition = "gad_total_score",mean = mean(tot_score), sd = sd(tot_score))

stats2 <- anx_gad_f2 %>% 
  semi_join(ms) %>% 
  summarise(definition = "gad_first_two",mean = mean(tot_score), sd = sd(tot_score))

gad_stats <- bind_rows(stats1,stats2)



