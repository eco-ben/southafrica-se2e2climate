library(StrathE2E2)
library(tidyverse)

plan(multisession, workers = 4)

model_path <- "../../StrathE2E_workspace/Models/"
variants <- list.dirs(str_glue("{model_path}South_Africa_MA/"), recursive = FALSE)
variants <- variants[!str_detect(variants, "2010-2015")] # Remove 2010-2015 variants
variants <- variants[str_detect(variants, "ssp")] # Include only variant folders with ESM/SSP info (remove 2010-2019 base variant)

variants <- str_split_i(variants, str_glue("{model_path}South_Africa_MA/"), 2)

load_and_run_variants <- function(model_variant, model_name, model_path, nyears = 50, output_file_base = "./outputs/initial_runs/final_biomass_", initial_cond_file_base = "./outputs/initial_runs/initial_conditions_") {
    print(str_glue("Running model for South Africa MA {model_variant}"))

    model <- e2e_read(model.name = model_name, model.variant = model_variant, models.path = model_path)
    results <- e2e_run(model, nyears = nyears, csv.output = FALSE)

    final_year <- StrathE2E2:::elt(results, "final.year.outputs")
    mass_results <- StrathE2E2:::elt(final_year, "mass_results_wholedomain")
    arrow::write_parquet(mass_results, paste0(output_file_base, model_variant, ".parq"))

    initial_results <- e2e_extract_start(model, results, csv.output = FALSE)
    initial_results$variable <- rownames(initial_results)
    names(initial_results) <- c("value", "variable")
    arrow::write_parquet(initial_results, paste0(initial_cond_file_base, model_variant, ".parq"))

    return()
}

walk(
    variants,
    function(x) {
        load_and_run_variants(x, "South_Africa_MA", model_path)
    }
)
