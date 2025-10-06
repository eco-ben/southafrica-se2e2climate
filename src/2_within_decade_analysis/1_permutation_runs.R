library(StrathE2E2)
library(tidyverse)
library(furrr)
library(glue)

plan(multisession, workers = 6)

models_path <- "../../StrathE2E_workspace/Models/"
output_path <- "./outputs/decade"

model <- e2e_read("South_Africa_MA", "2010-2019-CNRM-ssp126", models.path = models_path)
setup_files <- StrathE2E2:::pkg.env$SETUPFILES
decades <- c("2010-2019", "2030-2039", "2060-2069")

for (decade in decades) {
    variants_to_permute <- c(glue("{decade}-CNRM-ssp126"), glue("{decade}-CNRM-ssp370"), glue("{decade}-GFDL-ssp126"), glue("{decade}-GFDL-ssp370"))

    permutation_plan <- e2e_driver_permutation_plan(model, variants_to_permute)
    variable_group_levels <- permutation_plan$variable_group_unique_levels

    permutation_plan <- permutation_plan$plan
    write.csv(permutation_plan, paste0(output_path, decade, "_permutation_plan.csv"), row.names = FALSE)

    permutation_sources <- split(permutation_plan, permutation_plan$perm_id)

    future_map(
        permutation_sources,
        function(x) {
            # Due to StrathE2E2 architecture, the setup files are not copied across. We need to reassign them on the workers.
            assign("SETUPFILES", setup_files, envir = getNamespace("StrathE2E2")$pkg.env)
            e2e_run_driver_permutation(x, model)
        },
        .options = furrr_options(
            globals  = list(model = model, setup_files = setup_files),
            packages = "StrathE2E2"
        )
    )
}
