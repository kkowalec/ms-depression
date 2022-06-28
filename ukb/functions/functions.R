library(tidyverse)
library(data.table)
library(glue)



check_for_code <- function(codes, col_name, data){
  # A function for dealing with datasets where you want to find all individuals with a specific code,
  # where the code can be in any of the columns.
  # The function assumes that the first column is the identifier for each row
  
  
  derive <- function(data, code){
    if(code %>% is.character()){
      column <- data %>% 
        filter(if_any(-1, ~(str_detect(.,code)))) %>% 
        mutate(x = 1L) %>% select(1, x) 
      names(column)[2] <- code
      
    } else {
      column <- data %>% 
        filter(if_any(2:ncol(data),  ~.x %in% all_of(code))) %>% 
        mutate(x = 1L) %>% select(1, x) 
      names(column)[2] <- paste0("code", code)
    }
    
    print( paste0("For code: ", code, " found a total of ", nrow(column), " cases"))
    column %>% 
      mutate(temp = 1L) %>% 
      select(1, temp) %>% 
      tibble()
  }
  
  # Check if the input is one or multiple codes
  if(length(codes) == 1){
    derive(data, codes) %>% 
      mutate(temp = 1L) %>% 
      select(1, {{ col_name }} := temp) %>% 
      tibble()
    
  } else {
    # Sometimes you want to merge multiple codes into "one" definition.
    output <- 
      codes %>% 
      map(derive, data = data) %>% 
      reduce(full_join, by = "f.eid") %>% 
      mutate(temp = 1L) %>% 
      select(1, {{ col_name }} := temp) %>% 
      tibble()
    print(
      paste0("All joined together, found a total of ", nrow(output), " cases")
    )
    return(output)
  }
}



derive_ukb_n <- function(icd_codes = NULL, selfrep_codes = NULL, medicine_codes = NULL, merge = TRUE, col_name){
  
  
  stopifnot(all(exists(c("icd_data", "selfrep_data", "meds_data"))))
  
  output <- list()
  if(!is.null(icd_codes)){
    icd_phenotypes <- check_for_code(codes = icd_codes, col_name = col_name, data =  icd_data)
    names(icd_phenotypes)[[2]] <- paste0(names(icd_phenotypes)[[2]], "_icd")
    output <- append(output, icd_phenotypes)
  }
  
  if(!is.null(selfrep_codes)){
    selfrep_phenotypes <- 
      check_for_code(codes = selfrep_codes, col_name = col_name, data =  selfrep_data) %>% tibble()
    names(selfrep_phenotypes)[[2]] <- paste0(names(selfrep_phenotypes)[[2]], "_sr")
    
    output <- append(output, selfrep_phenotypes)
  }
  
  if(!is.null(medicine_codes)){
    selfrep_med <- check_for_code(codes = medicine_codes, col_name = col_name, data =  meds_data)
    names(selfrep_med)[[2]] <- paste0(names(selfrep_med)[[2]], "_med")
    output <- output %>% append(., selfrep_med)
  }
  
  list(icd_phenotypes, selfrep_phenotypes, selfrep_med)  %>% 
    reduce(full_join, by = "f.eid")
  
  
}



predict_with_prs <- function(dependant, predictor,data, binary=FALSE){
  if(binary){
    data <- data %>% 
      filter(.data[[dependant]] == 0 | .data[[dependant]] == 1)
  }
  
  predictors <- data %>% 
    select(all_of(predictor), starts_with("PC"))
  
  dependant_var <- data %>% 
    select(dependant)
  
  if (binary) {
    predictors %>% 
      glm(dependant_var[[1]] ~ ., data = .) %>% 
      tidy() %>% 
      filter(term %in% predictor) %>% 
      add_column(phenotype = dependant)
    
  } else {
    predictors %>% 
      lm(dependant_var[[1]] ~ ., data = .) %>% 
      tidy() %>% 
      filter(term %in% predictor) %>% 
      add_column(phenotype = dependant)
    
  }
}

probable_case <- function(tbl){
  tbl %>% 
    rowwise() %>%
    mutate(total = sum(c_across(2:4), na.rm = T)) %>% 
    ungroup() %>% 
    filter(total >= 2) %>% 
    inner_join(., qc_cohort, by = "f.eid")
  
}

####
# FUNCTION FOR EXTRACTING DATA FROM UKB ON VECTOR
check_dataset <- function(variable, dataset){
  ## Given a dataset filepath and a colname string, check if the variable exists as a column,
  ## and if it exists, fread it.
  
  # Can only do for one code at a time.
  stopifnot(length(variable) == 1)
  
  # Read in column names
  dataset_columns <- fread(dataset, nrows = 0) %>% colnames()
  
  # extract overlap
  cols <- dataset_columns[str_detect(dataset_columns,variable)]
  
  if(!rlang::is_empty(cols)){
    
    output <- fread(dataset, select = append("f.eid", cols))
    print(glue(
      "Found a total of {ncol(output)-1} columns using {variable} from:
      {dataset}"
    ))
    output
  } else {
    print(glue(
      "Found no columns in 
      {dataset}
      "))
    return(NULL)
  }
  
}

check_ukb <- function(variable, ukb_data){
  # wrapper for looping over several datasets
  # ukb_data is a hard coded filepaths for UKB data on vector
  # Filepaths for the UKB data
  ukb_data <- list("/nfs/AGE/UKB/data/211014/r/ukb48847.tab",
                   "/nfs/AGE/UKB/data/211019/r/ukb48978.tab",
                   "/nfs/AGE/UKB/data/180524/r/ukb22140.tab")
  
  map(ukb_data, variable = {{ variable }}, check_dataset) %>% 
    # remove empty entries
    compact() %>% 
    reduce(full_join, by = "f.eid")
  
}


extract_ukb <- function(variables){
  # if inputting multiple variables, iterate over each.
  if(length(variables) == 1){
    check_ukb( {{ variables }} ) %>% tibble()
    
  } else{
    map( {{ variables }} , check_ukb) %>% 
      reduce(full_join, by = "f.eid") %>% tibble()
  }
  
}