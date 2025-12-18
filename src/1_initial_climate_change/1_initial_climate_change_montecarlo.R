library(StrathE2E2)
library(tidyverse)
library(furrr)

n_workers <- 20
plan(multisession, workers = n_workers)

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

baseline_variants <- c("CNRM" = "2010-2019-CNRM-ssp126", "GFDL" = "2010-2019-GFDL-ssp126")
n_iter = 1000
base_result_directory <- "./outputs/initial_runs/monte_carlo_results/"

for (variant in c(variants[variants %in% baseline_variants], variants[!variants %in% baseline_variants])) {
    esm <- str_extract(variant, "[:upper:]{4}")
    base_parms_ident <- paste0(baseline_variants[esm], "-MC")

    baseline_mode <- FALSE
    if (variant %in% baseline_variants) {
        baseline_mode <- TRUE
    }

    variant_result_dir <- file.path(base_result_directory, "South_Africa_MA", variant, "CredInt")
    
    if (!dir.exists(variant_result_dir)) {
        dir.create(variant_result_dir, recursive=TRUE)
    }
    
    if (!baseline_mode) {
      baseline_param_files <- list.files(file.path(
        base_result_directory, 
        "South_Africa_MA", 
        baseline_variants[esm], 
        "CredInt"
      ), pattern = "parameters", full.names = TRUE)
      file.copy(baseline_param_files, variant_result_dir)
    }
    

    model_setup <- data.frame(
        region_name = "South_Africa_MA",
        model_variant = variant,
        model_path = model_path,
        result_path = base_result_directory
    )
    batch_setup <- do.call("rbind", replicate(n_workers, model_setup, simplify = FALSE))
    batch_setup$batch_id <- seq_len(n_workers)

    if (baseline_mode) {
        model_ident_label <- "baseline"
    } else {
        model_ident_label <- "scenario"
    }
    batch_setup$batch_model_ident <- paste(batch_setup$model_variant, "batch", batch_setup$batch_id, sep = "-")

    batch_setup$n_iter <- n_iter / n_workers
    batch_setup$begin_sample <- NA

    if (baseline_mode) {
        batch_setup$n_iter <- ifelse(batch_setup$batch_id == 1, batch_setup$n_iter, batch_setup$n_iter + 1)
    } else {
        for (r in seq_len(n_workers)) {
            if (r == 1) {
                batch_setup[r, ]$begin_sample <- 1
            } else {
                batch_setup[r, ]$begin_sample <- (batch_setup[r, ]$n_iter * batch_setup[r - 1, ]$batch_id + 1)
            }
        }
    }

    batch_setup$baseline_mode <- baseline_mode
    batch_setup$base_parms_ident <- base_parms_ident

    batch_list <- split(batch_setup, batch_setup$batch_id)
    
    print(paste0("Running Monte Carlo simulations for ", variant))
    
    future_map(
        batch_list,
        function(x) {
          baseline_mode <- x$baseline_mode
            model <- e2e_read(
                x$region_name,
                x$model_variant,
                models.path = x$model_path,
                results.path = x$result_path,
                model.ident = x$batch_model_ident
            )

            if (baseline_mode) {
                e2e_run_mc(
                    model,
                    baseline.mode = baseline_mode,
                    n_iter = x$n_iter,
                    csv.output = TRUE,
                    runtime.plot = FALSE,
                    postprocess = FALSE
                )
            } else {
                e2e_run_mc(
                    model,
                    baseline.mode = baseline_mode,
                    baseparms.ident = x$base_parms_ident,
                    begin.sample = x$begin_sample,
                    n_iter = x$n_iter,
                    csv.output = TRUE,
                    runtime.plot = FALSE,
                    postprocess = FALSE
                )
            }

            return()
        },
        .progress = TRUE,
        .options = furrr_options(seed = TRUE)
    )

    combined_model <- e2e_read(
        "South_Africa_MA",
        variant,
        models.path = model_path,
        results.path = base_result_directory,
        model.ident = paste0(variant, "-MC")
    )
    combined_data <- e2e_merge_sens_mc(
        combined_model,
        selection = "MC",
        ident.list = batch_setup$batch_model_ident,
        postprocess = FALSE,
        csv.output = TRUE
    )
    processed_data <- e2e_process_sens_mc(
        combined_model,
        selection = "MC",
        csv.output = TRUE
    )
    
    batch_files <- list.files(variant_result_dir, full.names = TRUE, pattern = "batch")
    file.remove(batch_files)
}
