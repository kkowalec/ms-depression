library(data.table)
library(fs)
library(glue)
library(tidyverse)

source("functions/functions.R")
source("functions/repro_functions.R")

meds_coding <- read_tsv("assets/coding4.tsv")



icd_list <- list(c("E780", "E782", "E784", "E785"),
                 c("I10", "I11", "I12", "I13", "I15"),
                 c("I20", "I21", "I22", "I23", "I24", "I25"),
                 c("I70", "I738", "I739"),
                 c("J40", "J41","J42", "J43", "J45", "J46"),
                 c("E10", "E11"),
                 c("C"),
                 c("G43"),
                 c("E02", "E03", "E05", "E06"),
                 c("M32"),
                 c("M15", "M16", "M17", "M18", "M19", "M47"),
                 c("M80", "M81", "M82"),
                 c("M79.7"),
                 c("N1"),
                 c("K25", "K21"),
                 c("K7"),
                 c("K58"),
                 c("G40", "G41"),
                 c("F32", "F33", "F34"),
                 c("F40", "F41"),
                 c("F30","F31"),
                 c("F20"),
                 c("F10", "F11", "F12", "F13", "F14", "F15", "F16", "F18", "F19", "K70"),
                 c("F50"),
                 c("F60"),
                 c("F90"),
                 c("K50", "K51", "G35", "M05", "M06", "K90", "L40", "M35","D51","G70", "M45", "M315", "M073"),
                 c("E28"),
                 c("S02", "S06", "S09"))


selfrep_list <- list(c("1473"), # high cholesterol
                     c("1065","1072"), # hypertension
                     c("1066", "1074","1075","1076","1077","1078","1079","1080"), # heart trouble
                     c("1067", "1087", "1088"), # disease of arteries 
                     c("1111", "1472", "1113", "1112"), # lung trouble
                     c("1220","1222", "1223"), # diabetes
                     c("x"),             # This is cancer, should be NA
                     c("1265"), # migraine
                     c("1226", "1225", "1428", "1522","1224"), # thyroid disease
                     c("1381"), # lupus
                     c("1465"), # osteorthritis
                     c("1309"), # osteoporosis
                     c("1542"), # fibromyalgia
                     c("1192","1193","1194","1405"), # kidney disease
                     c("1400", "1138", "1142"), # open sore or ulcer
                     c("1158", "1579", "1580"), # liver problems
                     c("1154"), # IBS
                     c("1264"), # epilepsy
                     c("1286"), # depression
                     c("1287"), # anxiety disorder
                     c("1291"), # bipolar disorder
                     c("1289"), # schizophrenia
                     c("1408", "1409", "1410", "1604"), # substance abuse
                     c("1470"), # eating disorder
                     c("x"), # personality disorders
                     c("x"), # adhd 
                     c("1461", "1261", "1464", "1456", "1453", "1382", "1462", "1463", "1331","1260", "1437","1313", "1377"),
                     c("1350"), # pcos 
                     c("1266")) # brain injury


phenotype_names <- c("high_cholesterol", 
                     "hypertension", 
                     "heart_trouble",
                     "disease_of_arteries",
                     "lung_trouble",
                     "diabetes_mellitus",
                     "cancer",
                     "migraine",
                     "thyroid_disease",
                     "lupus",
                     "osteoarthritis",
                     "osteoporosis",
                     "fibromyalgia",
                     "kidney_disease",
                     "open_sore_or_ulcer",
                     "liver_problems",
                     "ibs",
                     "epilepsy",
                     "depression",
                     "anxiety_disorder",
                     "bipolar_disorder",
                     "schizophrenia",
                     "suds",
                     "eating_disorder",
                     "personality_disorder",
                     "adhd",
                     "autoimmune_disorder",
                     "pcos",
                     "brain_injury")


joined_definitions <- tibble(icd = icd_list, phenotype = phenotype_names, selfrep = selfrep_list)
save(joined_definitions, file = "produced_data/exclusion_criteria_for_controls.RData")


# List of meds for each autoimmune disorder
# The list of medications used in multiple sclerosis, rheumatoid arthritis and inflammatory bowel disorder were # # # # derived from the Glansville et al (2021). (https://www.bpsgos.org/article/S2667-1743(21)00004-5/fulltext).




ms_meds <- c("interferon beta",
             "interferon beta-1a", 
             "interferon beta-1b",
             "interferon beta-1b product", 
             "avonex 6million iu injection (pdr for recon)+solvent",
             "betaferon 300micrograms injection (pdr for recon)+diluent", 
             "glatiramer",
             "copaxone 20mg injection (pdr for recon)")

ibd_meds <- c("prednisolone",
              "methotrexate",
              "mesalazine",
              "asacol 400mg e/c tablet",
              "pentasa sr 250mg m/r tablet",
              "asacol mr 400mg e/c tablet", 
              "prednisolone product",
              "balsalazide disodium", 
              "humira 40mg injection solution 0.8ml prefilled syringe",
              "olsalazine","adalimumab",
              "mtx - methotrexate", 
              "methylprednisolone", 
              "azt - azathioprine",
              "5asa - mesalazine")



ra_meds <- c("prednisolone",
             "methotrexate",
             "azathioprine", 
             "hydroxychloroquine", 
             "leflunomide",
             "prednisolone product",
             "humira 40mg injection solution 0.8ml prefilled syringe",
             "adalimumab",
             "mtx - methotrexate",
             "methylprednisolone",
             "azt - azathioprine", 
             "arava 10mg tablet", 
             "arava 20mg tablet", 
             "gold product", 
             "sodium aurothiomalate", 
             "gold sodium thiomalate")


ms_meds_codes <- 
  meds_coding %>% 
  filter(meaning %in% ms_meds) %>% .$coding


ibd_meds_codes <- 
  meds_coding %>% 
  filter(meaning %in% ibd_meds) %>% .$coding


ra_meds_codes <- 
  meds_coding %>% 
  filter(meaning %in% ra_meds) %>% .$coding


auto_immune_meds <- list(ms = ms_meds_codes, ibd = ibd_meds_codes, ra = ra_meds_codes)
save(auto_immune_meds, file = "produced_data/immune_disorder_meds.RData")