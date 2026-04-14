library(stringr)
library(StrathE2E2)
setwd("C:/users/grier/Documents/Projects/southafrica-se2e2climate/")

variants <- list.dirs("../../southafrica_paper/S_Benguela_MA/", recursive=FALSE, full.names = FALSE)
variants <- variants[!str_detect(variants, "2010-2015")]
old_out_dir <- "./outputs/initial_runs/monte_carlo_results/S_Benguela_MA/"
new_out_dir <- "./outputs/initial_runs/monte_carlo_results/S_Benguela_MA/"

old_rel_path <- "C:/users/grier/OneDrive - University of Strathclyde/Documents/repos/southafrica-se2e2climate/"

for (var in variants) {
  old_out <- file.path(old_rel_path, file.path(old_out_dir, var))
  files <- list.files(old_out, recursive = TRUE, full.names = TRUE)
  files_to_copy <- files[str_detect(files, pattern = "parameters") | str_detect(files, pattern = "mass") | str_detect(files, pattern = "AAM")]
  
  new_out <- file.path(new_out_dir, var)
  new_fns <- gsub(pattern = old_out, replacement = new_out, x = files_to_copy)
  
  if (!dir.exists(dirname(new_fns[1]))) {
    dir.create(dirname(new_fns[1]), recursive = TRUE)
  }
  
  file.copy(from = files_to_copy, to = new_fns, overwrite = TRUE)
}


setwd("C:/users/grier/OneDrive - University of Strathclyde/Documents/repos/southafrica-se2e2climate/")
out_dir <- "./outputs/initial_runs/monte_carlo_results/"

for (variant in variants) {
  model <- e2e_read("S_Benguela_MA", variant, models.path = "../../../../Documents/southafrica_paper/", results.path = out_dir, model.ident =str_glue("{variant}-MC"))
  cum_files <- list.files(file.path(out_dir, "S_Benguela_MA", variant, "CredInt"), full.names = TRUE, pattern = "cumulative")
  
  inshore_aam_store <- cum_files[str_detect(cum_files, "cumulative_inshoreaamass")]
  inshore_aam_store <- read.csv(inshore_aam_store)
  
  offshore_aam_store <- cum_files[str_detect(cum_files, "cumulative_offshoreaamass")]
  offshore_aam_store <- read.csv(offshore_aam_store)
  
  whole_aam_store <- cum_files[str_detect(cum_files, "cumulative_wholeaamass")]
  whole_aam_store <- read.csv(whole_aam_store)
  
  StrathE2E2:::CredInt_make_aamass_results(model, inshore_aam_store, offshore_aam_store, whole_aam_store, csv.output=TRUE, creds = c(0.025, 0.25, 0.5, 0.75, 0.975))
}
