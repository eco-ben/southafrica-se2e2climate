library(furrr)
library(data.table)

plan(multisession, workers = 4)
source("./src/model_permutations_common.R")

#### Create climate permutation sets for across decade analyses ----
within_esm_ssp_permutations <- read.csv("./outputs/within_esm_ssp_decade_permutations.csv")
within_esm_ssp_permutations$constants <- "2010-2019-CNRM-ssp126"

master <- load_all_possible_drivers("../../StrathE2E_workspace/Models/South_Africa_MA/", "South_Africa_MA")
unique_variables <- unique(master$variable)

within_esm_ssp_rows_list <- apply(within_esm_ssp_permutations, 1, function(row) as.list(row))
within_esm_ssp_permutation_sources <- future_map(within_esm_ssp_rows_list, create_variable_source)

within_esm_ssp_permutation_sources <- rbindlist(within_esm_ssp_permutation_sources)
arrow::write_parquet(within_esm_ssp_permutation_sources, "./outputs/within_esm_ssp_decade_permutation_sources.parq")
