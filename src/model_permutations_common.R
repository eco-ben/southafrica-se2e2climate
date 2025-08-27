source("src/build_model_forcing_functions.R")

master <- load_all_possible_drivers("../../StrathE2E_workspace/Models/South_Africa_MA/", "South_Africa_MA")

decades <- c("2010-2019", "2030-2039", "2060-2069") # Define decades that I want to investigate (currently only 2010-2019, 2030-39 and 2060-69)
ESMs <- unique(master$ESM)
SSPs <- unique(master$SSP)
blocks <- c(
    "light", "temperature", "river_outputs", "vertical_mixing", "water_flows", "nutrient_conc", "atm_nut_flux"
)

# Function to map numeric index to actual ESM/SSP pair
map_variant_index <- function(options, index) {
    paste(options[index, ], collapse = "-")
}

variable_groups <- list(
    constants = c(
        "so_logespm", "si_logespm", "Inshore_waveheight", "s1_pdist",
        "s2_pdist", "s3_pdist", "d1_pdist", "d2_pdist", "d3_pdist", "DO_mixlscale",
        "mixlscale", "si_othernitrate", "si_otherammonia", "so_othernitrate",
        "so_otherammonia", "DO_logkvert", "DO_d_upwelling", "d_DO_downwelling",
        "DO_nitrate", "DO_ammonia", "DO_detritus"
    ),
    light = "sslight",
    temperature = c("so_temp", "d_temp", "si_temp"),
    river_outputs = c("rivervol", "rivnitrate", "rivammonia", "rivdetritus"),
    vertical_mixing = c("logkvert", "d_so_upwelling", "so_d_downwelling"),
    water_flows = c("so_inflow", "d_inflow", "si_inflow", "si_outflow", "so_si_flow"),
    nutrient_conc = c(
        "so_nitrate", "so_ammonia", "so_phyt", "so_detritus",
        "d_nitrate", "d_ammonia", "d_phyt", "d_detritus",
        "si_nitrate", "si_ammonia", "si_phyt", "si_detritus"
    ),
    atm_nut_flux = c("so_atmnitrate", "so_atmammonia", "si_atmnitrate", "si_atmammonia")
)

get_variable_group <- function(var_name, var_groups = variable_groups) {
    for (group in names(var_groups)) {
        if (var_name %in% var_groups[[group]]) {
            return(group)
        }
    }
    return(NA) # If not found in any group
}

create_variable_source <- function(run_perm) {
    # Create empty variable source dataframe
    variable_sources <- data.frame(
        perm_id = run_perm["perm_id"],
        variable = unique_variables,
        group = sapply(unique_variables, get_variable_group),
        SSP = rep(NA, length(unique_variables)),
        ESM = rep(NA, length(unique_variables)),
        decade = rep(NA, length(unique_variables))
    )

    variable_sources$SSP <- sapply(variable_sources$group, function(x) get_var_info(x, run_perm, "[:lower:]{3}\\d{3}"))
    variable_sources$ESM <- sapply(variable_sources$group, function(x) get_var_info(x, run_perm, "[:upper:]{4}"))
    variable_sources$decade <- sapply(variable_sources$group, function(x) get_var_info(x, run_perm, "\\d{4}-\\d{4}"))

    return(variable_sources)
}

get_var_info <- function(v_group, run_perm, pattern) {
    info <- as.character(run_perm[v_group])
    return(str_extract(info, pattern))
}

save_final_year_output <- function(model_results) {
    final_year <- StrathE2E2:::elt(model_results, "final.year.outputs")
    mass_results <- StrathE2E2:::elt(final_year, "mass_results_wholedomain")
    return(mass_results)
}

run_permutation_model <- function(perm_model, perm_id, nyears = 50, output_summarisation, output_file_base = "./outputs/model_output_", initial_cond_file_base = "./outputs/initial_conditions_") {
    results <- e2e_run(perm_model, nyears = nyears, csv.output = FALSE)
    result_to_save <- output_summarisation(results)

    arrow::write_parquet(result_to_save, paste0(output_file_base, "perm_", perm_id, ".parq"))

    initial_results <- e2e_extract_start(perm_model, results, csv.output = FALSE)
    initial_results$variable <- rownames(initial_results)
    names(initial_results) <- c("value", "variable")
    arrow::write_parquet(initial_results, paste0(initial_cond_file_base, "perm_", perm_id, ".parq"))

    return(NULL)
}

load_and_run_perm <- function(master_variables, variable_sources, perm_id, nyears = 50, output_summarisation, output_file_base = "./outputs/model_output_", initial_cond_file_base = "./outputs/initial_conditions_", t_model_name = "South_Africa_MA", t_model_variant = "2010-2019-CNRM-ssp126", models_path = "../../StrathE2E_workspace/Models/") {
    t_model <- e2e_read(model.name = t_model_name, model.variant = t_model_variant, models.path = models_path)
    perm_model <- rebuild_model_drivers(t_model, master_variables, variable_sources)
    run_permutation_model(perm_model, perm_id, nyears = nyears, output_summarisation = output_summarisation, output_file_base = output_file_base, initial_cond_file_base = initial_cond_file_base)

    return(NULL)
}
