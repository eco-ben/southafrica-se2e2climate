library(StrathE2E2)
library(tidyverse)
library(furrr)
library(glue)

plan(multisession, workers = 20)

models_path <- "../../StrathE2E_workspace/Models/"
output_path <- "./outputs/across_esm_permutations/"

decades <- c("2010-2019", "2020-2029", "2030-2039", "2040-2049", "2050-2059", "2060-2069")
ssps <- c("ssp126", "ssp370")
esms <- c("CNRM", "GFDL")
model_parameterisations <- esms

for (model_parm in model_parameterisations) {
    model <- e2e_read("South_Africa_MA", glue("2010-2019-{model_parm}-ssp126"), models.path = models_path)
    setup_files <- StrathE2E2:::pkg.env$SETUPFILES

    for (decade in decades) {
        for (ssp in ssps) {
            sub_output_path <- file.path(output_path, glue("{model_parm}_parameterisation"), ssp, decade)
            variants_to_permute <- c(glue("{decade}-{esms[1]}-{ssp}"), glue("{decade}-{esms[2]}-{ssp}"))

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

            cat("Running permutations for", model_parm, ssp, decade)
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
                )
            )
        }
    }
}
