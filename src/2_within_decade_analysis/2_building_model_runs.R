library(furrr)
library(data.table)

plan(multisession, workers = 4)
source("./src/model_permutations_common.R")

#### Create climate permutation sets for within decade analyses ----
within_decade_permutations <- read.csv("./outputs/within_decade_ssp_esm_permutations.csv")
within_decade_permutations$constants <- "2010-2019-CNRM-ssp126"

master <- read.csv("./outputs/master_forcings_South_Africa_MA.csv")
unique_variables <- unique(master$variable)

within_decade_rows_list <- apply(within_decade_permutations, 1, function(row) as.list(row))
within_decade_permutation_sources <- future_map(within_decade_rows_list, create_variable_source)

within_decade_permutation_sources <- rbindlist(within_decade_permutation_sources)
arrow::write_parquet(within_decade_permutation_sources, "./outputs/within_decade_ssp_esm_permutation_sources.parq")
