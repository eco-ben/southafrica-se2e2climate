# se2e_workspace <- "../../StrathE2E_workspace"
# analysis_workspace <- "C:/Users/kbb25108/OneDrive - University of Strathclyde/Documents/repos/southafrica-se2e2climate/"

# setwd(se2e_workspace)

# se2e_physics_names <- c(
#     "month", "sslight", "so_logespm", "si_logespm",
#     "so_temp", "d_temp", "si_temp", "rivervol", "logkvert",
#     "mixlscale", "d_so_upwelling", "so_d_downwelling", "so_inflow",
#     "d_inflow", "si_inflow", "si_outflow", "so_si_flow",
#     "s1_pdist", "s2_pdist", "s3_pdist", "d1_pdist", "d2_pdist",
#     "d3_pdist", "Inshore_waveheight", "DO_logkvert", "DO_mixlscale",
#     "DO_d_upwelling", "d_DO_downwelling"
# )
# se2e_chemistry_names <- c(
#     "month",
#     "so_nitrate", "so_ammonia", "so_phyt", "so_detritus",
#     "d_nitrate", "d_ammonia", "d_phyt", "d_detritus",
#     "si_nitrate", "si_ammonia", "si_phyt", "si_detritus",
#     "rivnitrate", "rivammonia", "rivdetritus",
#     "so_atmnitrate", "so_atmammonia",
#     "si_atmnitrate", "si_atmammonia",
#     "si_othernitrate", "si_otherammonia",
#     "so_othernitrate", "so_otherammonia",
#     "DO_nitrate", "DO_ammonia", "DO_detritus"
# )
# variable_groups <- list(
#     constants = c(
#         "so_logespm", "si_logespm", "Inshore_waveheight", "s1_pdist",
#         "s2_pdist", "s3_pdist", "d1_pdist", "d2_pdist", "d3_pdist", "DO_mixlscale",
#         "mixlscale", "si_othernitrate", "si_otherammonia", "so_othernitrate",
#         "so_otherammonia", "DO_logkvert", "DO_d_upwelling", "d_DO_downwelling",
#         "DO_nitrate", "DO_ammonia", "DO_detritus"
#     ),
#     light = "sslight",
#     temperature = c("so_temp", "d_temp", "si_temp"),
#     river_outputs = c("rivervol", "rivnitrate", "rivammonia", "rivdetritus"),
#     vertical_mixing = c("logkvert", "d_so_upwelling", "so_d_downwelling"),
#     water_flows = c("so_inflow", "d_inflow", "si_inflow", "si_outflow", "so_si_flow"),
#     nutrient_conc = c(
#         "so_nitrate", "so_ammonia", "so_phyt", "so_detritus",
#         "d_nitrate", "d_ammonia", "d_phyt", "d_detritus",
#         "si_nitrate", "si_ammonia", "si_phyt", "si_detritus"
#     ),
#     atm_nut_flux = c("so_atmnitrate", "so_atmammonia", "si_atmnitrate", "si_atmammonia")
# )

# # Function to map numeric index to actual ESM/SSP pair
# map_variant_index <- function(options, index) {
#     if (class(options) == "data.frame") {
#         return(paste(options[index, ], collapse = "-"))
#     }
#     if (class(options) == "character") {
#         return(paste(options[index], collapse = "-"))
#     }
# }

# load_all_possible_drivers <- function(model_base_path, model_name, physics_names = se2e_physics_names, chemistry_names = se2e_chemistry_names) {
#     #' Load all possible physics and chemistry driving values from model variants.
#     #'
#     #' @description This function creates a master dataframe containing all variants of
#     #' physics and chemistry drivers. (Note that some drivers do not vary across model variants).
#     #'
#     #' @param model_base_path character. Base model folder that contains variant subfolders.
#     #' @param model_name character. Model name contained in driving data file names.
#     #' @param physics_names character. Correct names to apply to loaded physics file (copied from StrathE2E2).
#     #' @param chemistry_names character. Correct names to apply to loaded chemistry file (copied from StrathE2E2).
#     #' @return Longform dataframe containing all possible variant driving data for chemistry and physics variables.
#     #'
#     #' @references StrathE2E2
#     #' @examples
#     #' load_all_possible_drivers("Models/South_Africa_MA/", "South_Africa_MA")

#     possible_variants <- dir(model_base_path)
#     possible_variants <- possible_variants[str_detect(possible_variants, "-[:upper:]{4}-[:lower:]{3}\\d{3}")]
#     decades <- unique(str_extract(possible_variants, "\\d{4}-\\d{4}"))
#     ESMs <- unique(str_extract(possible_variants, "[:upper:]{4}"))
#     SSPs <- unique(str_extract(possible_variants, "[:lower:]{3}\\d{3}"))

#     master_forcings <- expand.grid(
#         month = 1:12,
#         decade = decades,
#         SSP = SSPs,
#         ESM = ESMs
#     )

#     for (variant in possible_variants) {
#         ESM <- str_extract(variant, "[:upper:]{4}")
#         SSP <- str_extract(variant, "[:lower:]{3}\\d{3}")
#         decade <- str_extract(variant, "\\d{4}-\\d{4}")

#         physics <- read.csv(paste0(
#             model_base_path,
#             decade, "-", ESM, "-", SSP,
#             "/Driving/",
#             "physics_", toupper(model_name), "_",
#             decade, "-", ESM, "-", SSP,
#             ".csv"
#         ))
#         names(physics) <- physics_names
#         physics$decade <- decade
#         physics$SSP <- SSP
#         physics$ESM <- ESM

#         chemistry <- read.csv(paste0(
#             model_base_path,
#             decade, "-", ESM, "-", SSP,
#             "/Driving/",
#             "chemistry_", toupper(model_name), "_",
#             decade, "-", ESM, "-", SSP,
#             ".csv"
#         ))
#         names(chemistry) <- chemistry_names
#         chemistry$decade <- decade
#         chemistry$SSP <- SSP
#         chemistry$ESM <- ESM

#         if (which(variant == possible_variants) == 1) {
#             master_forcings <- left_join(master_forcings, physics, by = c("month", "decade", "SSP", "ESM"))
#             master_forcings <- left_join(master_forcings, chemistry, by = c("month", "decade", "SSP", "ESM"))
#         } else {
#             master_forcings <- rows_update(master_forcings, physics, by = c("month", "decade", "SSP", "ESM"))
#             master_forcings <- rows_update(master_forcings, chemistry, by = c("month", "decade", "SSP", "ESM"))
#         }
#     }

#     master_forcings <- pivot_longer(
#         master_forcings,
#         cols = names(master_forcings)[!names(master_forcings) %in% c("month", "SSP", "ESM", "decade")],
#         names_to = "variable",
#         values_to = "value"
#     )
# }

# get_variable_group <- function(var_name, var_groups) {
#     for (group in names(var_groups)) {
#         if (var_name %in% var_groups[[group]]) {
#             return(group)
#         }
#     }
#     return(NA) # If not found in any group
# }

# get_var_info <- function(v_group, run_perm, pattern) {
#     info <- as.character(run_perm[v_group])
#     return(str_extract(info, pattern))
# }

# create_variable_source <- function(run_perm, unique_variables, variable_groups) {
#     #' Takes a the variable group sources for a single permutation and translates them into sources for each variable within each group
#     #'
#     #' @description The `rebuild_model_drivers()` function requires the input sources (variant source)
#     #' for each variable individually rather than grouped into variable groups. We do this
#     #' but taking the variable group allocations and expanding them

#     # Create empty variable source dataframe
#     variable_sources <- data.frame(
#         perm_id = run_perm["perm_id"],
#         variable = unique_variables,
#         group = sapply(unique_variables, function(x) get_variable_group(x, variable_groups)),
#         SSP = rep(NA, length(unique_variables)),
#         ESM = rep(NA, length(unique_variables)),
#         decade = rep(NA, length(unique_variables))
#     )

#     variable_sources$SSP <- sapply(variable_sources$group, function(x) get_var_info(x, run_perm, "[:lower:]{3}\\d{3}"))
#     variable_sources$ESM <- sapply(variable_sources$group, function(x) get_var_info(x, run_perm, "[:upper:]{4}"))
#     variable_sources$decade <- sapply(variable_sources$group, function(x) get_var_info(x, run_perm, "\\d{4}-\\d{4}"))

#     return(variable_sources)
# }

# rebuild_model_drivers <- function(template_model, master_forcing_values, variable_sources, physics_names = se2e_physics_names, chemistry_names = se2e_chemistry_names) {
#     #' Rebuild model drivers from different sources.
#     #'
#     #' @description Rebuild template_model drivers using different variant values from variable_sources
#     #'
#     #' @param template_model StrathE2E2 model. Template e2e_read() model to recreate model structure.
#     #' @param master_forcing_values dataframe. Master dataframe containing all possible variants of driving data.
#     #' (col-names = month, decade, SSP, ESM, variable, value).
#     #' @param variable_sources dataframe. Dataframe containing desired value sources for each driving variable.
#     #' (col-names = variable, SSP, ESM, decade).
#     #' @param physics_names character. Physics variable names (copied from StrathE2E2).
#     #' @param chemistry_names character. Chemistry variable names (copied from StrathE2E2).
#     #' @return StrathE2E2 model object containing driving data from variant sources indicated in variable_sources.
#     #'
#     #' @references StrathE2E2
#     #' @examples
#     #' rebuild_model_drivers(
#     #'     e2e_read(model.name = "South_Africa_MA", model.variant = "2010-2019-CNRM-ssp126"),
#     #'     load_all_possible_drivers("./Models/South_Africa_MA/", "South_Africa_MA"),
#     #'     variable_sources
#     #' )

#     data <- StrathE2E2:::elt(template_model, "data")
#     current_physics <- StrathE2E2:::elt(data, "physics.drivers")
#     current_chemistry <- StrathE2E2:::elt(data, "chemistry.drivers")

#     new_physics <- current_physics
#     for (p_name in physics_names) {
#         if (p_name == "month") {
#             next
#         }
#         var_source <- variable_sources[variable_sources$variable == p_name, ]

#         if (paste0(var_source$decade, "-", var_source$ESM, "-", var_source$SSP) == template_model$setup$model.variant) {
#             next
#         }

#         new_vals <- master_forcing_values[
#             master_forcing_values$decade == var_source$decade &
#                 master_forcing_values$ESM == var_source$ESM &
#                 master_forcing_values$SSP == var_source$SSP &
#                 master_forcing_values$variable == p_name,
#         ]

#         new_physics[, p_name] <- new_vals[sort(new_vals$month), ]$value
#     }

#     new_chemistry <- current_chemistry
#     for (c_name in chemistry_names) {
#         if (c_name == "month") {
#             next
#         }
#         var_source <- variable_sources[variable_sources$variable == c_name, ]

#         if (paste0(var_source$decade, "-", var_source$ESM, "-", var_source$SSP) == template_model$setup$model.variant) {
#             next
#         }

#         new_vals <- master_forcing_values[
#             master_forcing_values$decade == var_source$decade &
#                 master_forcing_values$ESM == var_source$ESM &
#                 master_forcing_values$SSP == var_source$SSP &
#                 master_forcing_values$variable == c_name,
#         ]
#         new_physics[, c_name] <- new_vals[sort(new_vals$month), ]$value
#     }

#     new_data <- data
#     new_data[["physics.drivers"]] <- new_physics
#     new_data[["chemistry.drivers"]] <- new_chemistry

#     new_model <- list(
#         setup = template_model$setup,
#         data = new_data
#     )

#     return(new_model)
# }

# save_final_year_output <- function(model_results) {
#     final_year <- StrathE2E2:::elt(model_results, "final.year.outputs")
#     mass_results <- StrathE2E2:::elt(final_year, "mass_results_wholedomain")
#     return(mass_results)
# }

# run_permutation_model <- function(perm_model, perm_id, nyears = 50, output_summarisation, output_file_base = "./outputs/model_output_", initial_cond_file_base = "./outputs/initial_conditions_") {
#     results <- e2e_run(perm_model, nyears = nyears, csv.output = FALSE)
#     result_to_save <- output_summarisation(results)

#     arrow::write_parquet(result_to_save, paste0(output_file_base, "perm_", perm_id, ".parq"))

#     initial_results <- e2e_extract_start(perm_model, results, csv.output = FALSE)
#     initial_results$variable <- rownames(initial_results)
#     names(initial_results) <- c("value", "variable")
#     arrow::write_parquet(initial_results, paste0(initial_cond_file_base, "perm_", perm_id, ".parq"))

#     return(NULL)
# }

# load_and_run_perm <- function(master_variables, variable_sources, perm_id, nyears = 50, output_summarisation, output_file_base = "./outputs/model_output_", initial_cond_file_base = "./outputs/initial_conditions_", t_model_name = "South_Africa_MA", t_model_variant = "2010-2019-CNRM-ssp126", models_path = "../../StrathE2E_workspace/Models/") {
#     #' Rebuild model drivers based on the driver sources in `variable_sources`, then run the altered model and save the required outputs and the initial conditions of the model.
#     #'
#     #' @param master_variables data.frame. Dataframe containing the required values to rebuild
#     #' the template model with new driving values.
#     #' @param variable_sources data.frame. Dataframe containing the new sources for driving
#     #' values that are going to be taken from `master_variables` and put into the template model.
#     #' @param perm_id integer. Permutation ID (usually an integer number e.g. permutation 1).
#     #' @param nyears integer. Number of years to run the model for, default is 50.
#     #' @param output_summarisation function. A function to apply to the StrathE2E model_results
#     #' object that will return a dataframe that is saved out as the desired model outputs.
#     #' @param output_file_base character. The base filename to write the output file to using
#     #' the results from `output_summarisation` function.
#     #' @param initial_cond_file_base character. The base filename to write the model initial
#     #' conditions at steady state.
#     #' @param t_model_name character. The model name required for loading the template model.
#     #' @param t_model_variant character. The model variant required for loading the template model.
#     #' @param models_path character. The models folder path required for loading the template model object.
#     #'
#     #' @returns Nothing (All outputs are saved as files only).

#     # Create a template model to rebuild the drivers from (this could be externalised so it doesn't happen in each call of load_and_run_perm(), but it is quite quick anyway).
#     t_model <- e2e_read(model.name = t_model_name, model.variant = t_model_variant, models.path = models_path)

#     # Rebuild the template model, replacing the necessary drivers with the new drivers from `variable_sources` and `master_variables`.
#     perm_model <- rebuild_model_drivers(t_model, master_variables, variable_sources)

#     # Run the permutation model, saving the necessary outputs and the initial conditions.
#     run_permutation_model(perm_model, perm_id, nyears = nyears, output_summarisation = output_summarisation, output_file_base = output_file_base, initial_cond_file_base = initial_cond_file_base)

#     return(NULL)
# }


# ### Workflow: example - South_Africa_MA climate decade permutation for a single ESM-SSP
# master <- load_all_possible_drivers("../../StrathE2E_workspace/Models/South_Africa_MA/", "South_Africa_MA")

# decades <- c("2010-2019", "2030-2039", "2060-2069") # Define decades that I want to investigate (currently only 2010-2019, 2030-39 and 2060-69)
# ESMs <- unique(master$ESM)
# SSPs <- unique(master$SSP)
# blocks <- c(
#     "light", "temperature", "river_outputs", "vertical_mixing", "water_flows", "nutrient_conc", "atm_nut_flux"
# )
# esm_ssp <- "CNRM-ssp126"

# # Create hybrid combinations within one ESM-SSP combo: 3^7 = 2,187 runs
# # 7 = number of variable groups that are being permuted
# single_esm_ssp_combinations <- expand.grid(
#     replicate(length(blocks), seq_len(length(decades)), simplify = FALSE),
#     stringsAsFactors = FALSE
# )

# single_esm_ssp_combinations <- as.data.frame(
#     t(apply(single_esm_ssp_combinations, 1, function(row) {
#         sapply(as.numeric(row), function(x) map_variant_index(decades, x))
#     }))
# )
# names(single_esm_ssp_combinations) <- blocks

# hybrid_sources <- data.frame(matrix(
#     NA,
#     nrow = nrow(single_esm_ssp_combinations),
#     ncol = length(blocks) + 1
# ))
# names(hybrid_sources) <- c("perm_id", names(single_esm_ssp_combinations))
# perm_id <- 1

# for (i in seq_len(nrow(single_esm_ssp_combinations))) {
#     row <- single_esm_ssp_combinations[i, ]
#     row <- sapply(row, function(x) paste0(x, "-", esm_ssp))

#     hybrid_sources[perm_id, ] <- c("perm_id" = perm_id, row)
#     perm_id <- perm_id + 1
# }

# write.csv(hybrid_sources, "./outputs/within_esm_ssp_decade_permutations.csv")
# # This output can be the output of the sampling function.

# # Now with these variant sources for each variable we need to translate those into sources
# # for each variable within the variable groups and pass that into the model rebuilding and
# # running functions

# within_esm_ssp_permutations <- hybrid_sources
# within_esm_ssp_permutations$constants <- "2010-2019-CNRM-ssp126"

# master <- load_all_possible_drivers("../../StrathE2E_workspace/Models/South_Africa_MA/", "South_Africa_MA")
# unique_variables <- unique(master$variable)

# within_esm_ssp_rows_list <- apply(within_esm_ssp_permutations, 1, function(row) as.list(row))
# within_esm_ssp_permutation_sources <- map(within_esm_ssp_rows_list, function(x) create_variable_source(x, unique_variables, variable_groups))
# within_esm_ssp_permutation_sources <- rbindlist(within_esm_ssp_permutation_sources)

# # Split the permutation sources dataframe into a list of source dataframes (one for each permutation run)
# # that can then be handled with the load and run functions.
# within_esm_ssp_permutation_sources <- split(within_esm_ssp_permutation_sources, within_esm_ssp_permutation_sources$perm_id)
# future_walk2(
#     within_esm_ssp_permutation_sources[1:200],
#     seq_len(length(within_esm_ssp_permutation_sources[1:200])),
#     function(x, y) {
#         load_and_run_perm(
#             master,
#             variable_sources = x,
#             perm_id = y,
#             output_summarisation = save_final_year_output,
#             output_file_base = "./outputs/across_decade_permutations/model_outputs_",
#             initial_cond_file_base = "./outputs/across_decade_permutations/init_conditions_"
#         )
#     }
# )

# # After model runs are complete then the outputs need to be collated and passed into the Shapley analysis function:
# permutations <- hybrid_sources # The design of all of the permutations
# result_files <- sapply(seq_len(nrow(hybrid_sources)), function(x) paste0("./outputs_across_decade_permutations/model_outputs_perm_", x, ".parq"))

# # For each permutation take the output file and extract the demersal fish biomass value
# results <- lapply(result_files, function(x) {
#     result <- read_parquet(x)
#     result <- result[result$Description == "Demersal fish", Model_annual_mean]
#     return(result)
# })
# results <- rbind(results)
# permutations$demersal_fish <- results

# input_variables <- names(permutations)[names(permutations) != "perm_id"]
# shapley_effects <- shapley_main_and_interactions(permutations, input_variables, "demersal_fish")
target_guilds <- c(
    "Total_nitrogen_mass",
    "Surface_layer_phytoplankton",
    "Deep_layer_phytoplankton",
    "Omnivorous_zooplankton",
    "Carnivorous_zooplankton",
    "Benthos_susp/dep_feeders",
    "Benthos_carn/scav_feeders",
    "Planktivorous_fish",
    "Migratory_fish",
    "Demersal_fish",
    "Birds",
    "Pinnipeds",
    "Cetaceans"
)
