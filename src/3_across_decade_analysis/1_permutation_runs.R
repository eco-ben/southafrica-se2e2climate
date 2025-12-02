library(StrathE2E2)
library(tidyverse)
library(furrr)
library(glue)

plan(multisession, workers = 20)

models_path <- "../../StrathE2E_workspace/Models/"
output_path <- "./outputs/across_decade_permutations/"

ESM_SSP <- c("CNRM-ssp126", "CNRM-ssp370", "GFDL-ssp126", "GFDL-ssp370")

for (esm_ssp in ESM_SSP) {
    model <- e2e_read("South_Africa_MA", glue("2010-2019-{esm_ssp}"), models.path = models_path)
    setup_files <- StrathE2E2:::pkg.env$SETUPFILES

    sub_output_path <- paste0(output_path, esm_ssp, "/")
    variants_to_permute <- c(glue("2010-2019-{esm_ssp}"), glue("2030-2039-{esm_ssp}"), glue("2060-2069-{esm_ssp}"))

    permutation_plan <- e2e_driver_permutation_plan(model, variants_to_permute)
    if (is.character(permutation_plan)) {
        for (variable in permutation_plan) {
            model$data$climate.matching[model$data$climate.matching$Driver == variable, ]$Grouping <- 0
        }
        permutation_plan <- e2e_driver_permutation_plan(model, variants_to_permute)
    }
    variable_group_levels <- permutation_plan$variable_group_unique_levels

    permutation_plan <- permutation_plan$plan
    write.csv(permutation_plan, paste0(sub_output_path, "permutation_plan.csv"), row.names = FALSE)

    permutation_sources <- split(permutation_plan, permutation_plan$perm_id)

    cat("Running permutations for", esm_ssp)
    future_map(
        permutation_sources,
        function(x) {
            # Due to StrathE2E2 architecture, the setup files are not copied across. We need to reassign them on the workers.
            assign("SETUPFILES", setup_files, envir = getNamespace("StrathE2E2")$pkg.env)
            e2e_run_driver_permutation(x, model, output_dir = sub_output_path)
        },
        .options = furrr_options(
            globals  = list(model = model, setup_files = setup_files, sub_output_path = sub_output_path),
            packages = "StrathE2E2"
        ),
        .progress = TRUE
    )
}

# Extracting ecosystem_indices using initial conditions

extract_ecosystem_indices <- function(model_results) {
    final_year <- StrathE2E2:::elt(model_results, "final.year.outputs")
    ecosystem_indices <- StrathE2E2:::elt(final_year, "ecosystem_indices")
    return(ecosystem_indices)
}

for (esm_ssp in ESM_SSP) {
    model <- e2e_read("South_Africa_MA", glue("2010-2019-{esm_ssp}"), models.path = models_path)
    setup_files <- StrathE2E2:::pkg.env$SETUPFILES

    sub_output_path <- paste0(output_path, esm_ssp, "/")
    permutation_plan <- read.csv(file.path(sub_output_path, "permutation_plan.csv"))
    names(permutation_plan) <- as.character(c(0:7, "perm_id"))

    permutation_sources <- split(permutation_plan, permutation_plan$perm_id)

    cat("Extracting ecosystem indices for", esm_ssp)
    future_map(
        permutation_sources,
        function(x) {
            # Attach initial conditions
            initial_conditions <- read.csv(file.path(sub_output_path, glue("init_conditions_perm_{unique(x$perm_id)}.csv")))
            model_init_cond <- model$data$initial.state
            for (i in seq_len(nrow(initial_conditions))) {
                model_init_cond[[i]] <- initial_conditions[i, ]$value
            }
            model$data$initial.state <- model_init_cond

            # Due to StrathE2E2 architecture, the setup files are not copied across. We need to reassign them on the workers.
            assign("SETUPFILES", setup_files, envir = getNamespace("StrathE2E2")$pkg.env)
            e2e_run_driver_permutation(
              x,
              model,
              nyears = 1,
              output_dir = sub_output_path,
              output_file_base = "ecosystem_indices_",
              output_summarisation = extract_ecosystem_indices
            )
        },
        .options = furrr_options(
          globals  = list(model = model, setup_files = setup_files, sub_output_path = sub_output_path, extract_ecosystem_indices=extract_ecosystem_indices),
          packages = c("StrathE2E2", "glue")
        ),
        .progress = TRUE
    )
}
