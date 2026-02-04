library(StrathE2E2)
library(tidyverse)
library(furrr)

n_workers <- 10
plan(multisession, workers = n_workers)

model_path <- "../../southafrica_paper/"
variants <- list.dirs(str_glue("{model_path}S_Benguela_MA/"), recursive = FALSE)
variants <- variants[!str_detect(variants, "2010-2015")] # Remove 2010-2015 variants
variants <- variants[str_detect(variants, "ssp")] # Include only variant folders with ESM/SSP info (remove 2010-2019 base variant)

variants <- str_split_i(variants, str_glue("{model_path}S_Benguela_MA/"), 2)
results_path_base <- "./outputs/initial_runs/"
sapply(variants, function(x) dir.create(file.path(results_path_base, x, "demfish_sens")))

params <- c("xxdfish", "xxdfishlar", "h_fishd", "h_fishdlar", "u_fishd", "u_fishdlar")
esms <- c("GFDL", "CNRM")
base_params <- sapply(
    esms,
    function(esm) {
        model_e <- e2e_read("S_Benguela_MA", str_glue("2010-2019-{esm}-ssp126"), models.path = model_path)
        sapply(params, function(x) {model_e$data$fitted.parameters[[x]]})
    }
)
base_params_save <- as.data.frame(base_params)
base_params_save$param_name <- rownames(base_params)
write.csv(base_params_save, "./outputs/initial_runs/demfish_sens_base_params.csv", row.names = FALSE)

variants <- c("2010-2019-CNRM-ssp370", "2010-2019-GFDL-ssp370") # Run only a subset of variants due to computational constraints
variation <- 0.25
steps <- 100

sampled_params <- lapply(
    esms,
    function(esm) {
        base_p <- base_params[, esm]
        lapply(params, function(x) {
            base_val <- base_p[x]
            sampled_p <- seq(base_val * (1 - variation), base_val * (1 + variation), length.out = steps)
            return(data.frame(param_name = rep(x, steps), param_value = sampled_p))
        })
    }
)
sampled_params <- lapply(sampled_params, function(x) {do.call(rbind, x)})
for(i in 1:length(sampled_params)) {sampled_params[[i]]$esm <- esms[i]}

param_runs <- do.call(rbind, sampled_params)
param_runs$param_id <- rep(seq(1, length(params) * steps), 2)
param_runs <- pivot_wider(param_runs, id_cols = c(param_name, param_id), names_from = esm, values_from = param_value)
write.csv(param_runs, file.path(results_path_base, "demfish_sens_params.csv"))

param_run <- function(variant, param_row) {
    esm <- str_extract(variant, "[:upper:]{4}")

    model_s <- e2e_read("S_Benguela_MA", variant, models.path = model_path)
    model_s$data$fitted.parameters[[param_row$param_name]] <- param_row[, esm]

    results <- e2e_run(model_s, nyears = 50, csv.output =  FALSE)

    initial <- e2e_extract_start(model_s, results, csv.output = FALSE)
    write.csv(initial, file.path(results_path_base, variant, "demfish_sens", str_glue("initial_conditions_dfishsens_p{param_row$param_id}.csv")))

    final_year <- results$final.year.outputs$mass_results_wholedomain
    write.csv(final_year, file.path(results_path_base, variant, "demfish_sens", str_glue("final_year_dfishsens_p{param_row$param_id}.csv")))
}

for (variant in variants) {
    param_list <- split(param_runs, param_runs$param_id)
    future_walk(param_list, function(x) {param_run(variant, x)}, .progress = TRUE)
}
