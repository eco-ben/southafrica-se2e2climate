library(furrr)
library(StrathE2E2)

plan(multisession, workers = 4)

source("./src/model_permutations_common.R")

within_esm_ssp_permutation_sources <- arrow::read_parquet("./outputs/within_esm_ssp_decade_permutation_sources.parq")
within_esm_ssp_permutation_sources <- split(within_esm_ssp_permutation_sources, within_esm_ssp_permutation_sources$perm_id)

future_walk2(
    within_esm_ssp_permutation_sources[1:200],
    seq_len(length(within_esm_ssp_permutation_sources[1:200])),
    function(x, y) {
        load_and_run_perm(
            master,
            variable_sources = x,
            perm_id = y,
            output_summarisation = save_final_year_output,
            output_file_base = "./outputs/across_decade_permutations/model_outputs_",
            initial_cond_file_base = "./outputs/across_decade_permutations/init_conditions_"
        )
    }
)
