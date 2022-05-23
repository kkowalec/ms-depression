save_results_on_git_hash <- function(folder, results){
  small_hash <- str_sub(system("git rev-parse HEAD", intern=TRUE)[1], start = -10)
  
  dir_path <- glue(folder,"/", small_hash)
  file <- Sys.Date() %>% str_replace_all(., "-", "_")
  filepath <- glue(dir_path, "/", file, ".RData")
  
  if(!dir_exists(dir_path)){
    print(glue("Folder {dir_path} does not exist. Creating it"))
    dir_create(dir_path)
  }
  
  filepath
  
  if(file_exists(filepath)){
    stop("A results file already exists on this commit")
  } 
  else{
    git_results <- list(results, params, sessionInfo())
    save(git_results ,file = filepath)
    print("Saved results, params and session info")
    
    
  }
  
}

get_results_on_git_hash <- function(folder){
  
  small_hash <- str_sub(system("git rev-parse HEAD", intern=TRUE)[1], start = -10)
  dir_path <- glue(folder,"/", small_hash)
  file <- Sys.Date() %>% str_replace_all(., "-", "_")
  glue(dir_path, "/", file, ".RData")
  
  
  
  
  
}