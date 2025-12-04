library(StrathE2E2)
library(tidyverse)
library(furrr)

plan(multisession, workers = 4)

model_path <- "../../StrathE2E_workspace/Models/"
variants <- list.dirs(str_glue("{model_path}South_Africa_MA/"), recursive = FALSE)
variants <- variants[!str_detect(variants, "2010-2015")] # Remove 2010-2015 variants
variants <- variants[str_detect(variants, "ssp")] # Include only variant folders with ESM/SSP info (remove 2010-2019 base variant)

variants <- str_split_i(variants, str_glue("{model_path}South_Africa_MA/"), 2)
results_path_base <- "./outputs/initial_runs/"
sapply(variants, function(x) dir.create(file.path(results_path_base, x)))


load_and_run_variants <- function(
    model_variant,
    model_name,
    model_path,
    nyears = 50,
    output_file_base = "final_biomass_",
    initial_cond_file_base = "initial_conditions_",
    indices_file_base = "ecosystem_indices_",
    flux_mat_file_base = "flux_matrix_exclspawn_",
    opt_res_file_base = "opt_results_",
    ann_flux_stats_file_base = "annual_flux_stats_",
    trophic_prod_file_base = "trophic_prod_") {
    print(str_glue("Running model for South Africa MA {model_variant}"))

    model <- e2e_read(model.name = model_name, model.variant = model_variant, models.path = model_path)
    results <- e2e_run(model, nyears = nyears, csv.output = FALSE)

    final_year <- StrathE2E2:::elt(results, "final.year.outputs")
    mass_results <- StrathE2E2:::elt(final_year, "mass_results_wholedomain")
    write.csv(mass_results, paste0(output_file_base, model_variant, ".csv"))

    network_results <- StrathE2E2:::elt(final_year, "ecosystem_indices")
    write.csv(network_results, paste0(indices_file_base, model_variant, ".csv"))

    flux_matrix <- StrathE2E2:::elt(final_year, "flow_matrix_excl_spawn_recruit")
    write.csv(flux_matrix, paste0(flux_mat_file_base, model_variant, ".csv"))

    opt_results <- StrathE2E2:::elt(final_year, "opt_results")
    write.csv(opt_results, paste0(opt_res_file_base, model_variant, ".csv"))

    annual_flux_stats <- StrathE2E2:::elt(final_year, "annual_flux_results_wholedomain")
    write.csv(annual_flux_stats, paste0(ann_flux_stats_file_base, model_variant, ".csv"))

    trophic_prod <- StrathE2E2:::elt(final_year, "HANPP_results")
    trophic_prod <- StrathE2E2:::elt(trophic_prod, "trophic_data")
    write.csv(trophic_prod, paste0(trophic_prod_file_base, model_variant, ".csv"))

    initial_results <- e2e_extract_start(model, results, csv.output = FALSE)
    initial_results$variable <- rownames(initial_results)
    names(initial_results) <- c("value", "variable")
    write.csv(initial_results, paste0(initial_cond_file_base, model_variant, ".csv"))

    return()
}

load_and_run_montecarlo <- function(model_variant, model_name, model_path, nyears = 50, n_iter = 1000, result_path = ".") {
    print(str_glue("Running monte carlo for South Africa MA {model_variant}"))

    model <- e2e_read(
        model.name = model_name,
        model.variant = model_variant,
        models.path = model_path,
        results.path = result_path
    )
    results <- e2e_run_mc(model, nyears = nyears, n_iter = n_iter, csv.output = TRUE)


    return()
}

future_map(
    variants,
    function(x) {
        load_and_run_variants(
            x,
            "South_Africa_MA",
            model_path,
            nyears = 1,
            output_file_base = file.path(results_path_base, x, "final_biomass_"),
            initial_cond_file_base = file.path(results_path_base, x, "initial_conditions_"),
            indices_file_base = file.path(results_path_base, x, "ecosystem_indices_"),
            flux_mat_file_base = file.path(results_path_base, x, "flux_matrix_exclspawn_"),
            opt_res_file_base = file.path(results_path_base, x, "opt_results_"),
            ann_flux_stats_file_base = file.path(results_path_base, x, "ann_flux_stats_"),
            trophic_prod_file_base = file.path(results_path_base, x, "trophic_prod_")
        )
    },
    .progress = TRUE
)

future_map(
    variants,
    function(x) {
        load_and_run_montecarlo(x, "South_Africa_MA", model_path, result_path = paste0(results_path_base, x))
    },
    .progress = TRUE
)
