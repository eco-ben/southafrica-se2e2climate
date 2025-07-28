source("./src/1_rebuilding_model_forcing.R")
library(furrr)

variable_groups <- list(
    constants = c("so_logespm", "si_logespm", "Inshore_waveheight", "s1_pdist", "s2_pdist", "s3_pdist", "d1_pdist", "d2_pdist", "d3_pdist"),
    light = "sslight",
    temperature = c("so_temp", "d_temp", "si_temp"),
    river_outputs = c("rivervol", "rivnitrate", "rivammonia", "rivdetritus"),
    vertical_mixing = c("logkvert", "mixlscale", "d_so_upwelling", "so_d_downwelling", "DO_logkvert", "DO_mixlscale", "DO_d_upwelling", "d_DO_downwelling"),
    water_flows = c("so_inflow", "d_inflow", "si_inflow", "si_outflow", "so_si_flow"),
    nutrient_conc = c(
        "so_nitrate", "so_ammonia", "so_phyt", "so_detritus",
        "d_nitrate", "d_ammonia", "d_phyt", "d_detritus",
        "si_nitrate", "si_ammonia", "si_phyt", "si_detritus"
    ),
    atm_nut_flux = c("so_atmnitrate", "so_atmammonia", "si_atmnitrate", "si_atmammonia"),
    other_nut_flux = c("si_othernitrate", "si_otherammonia", "so_othernitrate", "so_otherammonia", "DO_nitrate", "DO_ammonia", "DO_detritus")
)
get_variable_group <- function(var_name, var_groups = variable_groups) {
    for (group in names(var_groups)) {
        if (var_name %in% var_groups[[group]]) {
            return(group)
        }
    }
    return(NA) # If not found in any group
}

run_permutations <- read.csv("./outputs/climate_variant_permutations.csv")
run_permutations$constants <- "2010-2019-CNRM-ssp126"

master <- load_all_possible_drivers("../../StrathE2E_workspace/Models/South_Africa_MA/", "South_Africa_MA")

get_var_info <- function(v_group, run_perm, pattern) {
    info <- as.character(run_perm[v_group])
    return(str_extract(info, pattern))
}

unique_variables <- unique(master$variable)
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

plan(multisession, workers = 6)
run_rows_list <- apply(run_permutations, 1, function(row) as.list(row))
run_permutation_sources <- future_map(run_rows_list, create_variable_source)
