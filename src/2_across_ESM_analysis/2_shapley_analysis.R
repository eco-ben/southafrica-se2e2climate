library(StrathE2E2)
library(tidyverse)
library(glue)

output_path <- "./outputs/across_esm_permutations/"

decades <- c("2010-2019", "2030-2039", "2060-2069")
variable_group_names <- c(
    "constant",
    "light",
    "temperature",
    "river_outputs",
    "vertical_mixing",
    "water_flows",
    "nutrient_concentrations"
)
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

decades <- c("2010-2019", "2020-2029", "2030-2039", "2040-2049", "2050-2059", "2060-2069")
ssps <- c("ssp126", "ssp370")
esms <- c("CNRM", "GFDL")
model_parameterisations <- esms

# Add Shapley effects for net primary production
extract_ecosystem_indices_output <- function(csv_file, output_id) {
    result <- read.csv(csv_file)
    result <- result[result$Description == output_id, "Value"]

    return(result)
}

for (model_parm in model_parameterisations) {
    for (decade in decades) {
        for (ssp in ssps) {
            sub_output_path <- file.path(output_path, glue("{model_parm}_parameterisation"), ssp, decade)
            cat("Assessing Shapley Effects for", model_parm, ssp, decade)

            if (!file.exists(file.path(sub_output_path, "permutation_plan.csv"))) {
                file.copy(from = file.path(output_path, glue("{model_parm}_parameterisation"), ssp, glue("{decade}permutation_plan.csv")), to = file.path(sub_output_path, "permutation_plan.csv"))
                file.remove(file.path(output_path, glue("{model_parm}_parameterisation"), ssp, glue("{decade}permutation_plan.csv")))
            }

            shapley_results <- expand.grid(variable_group = variable_group_names, output = c(target_guilds, "netprimprod"), stringsAsFactors = FALSE)
            shapley_results$shapley_effect <- NA
            permutation_plan <- read.csv(file.path(sub_output_path, "permutation_plan.csv"), col.names = c(variable_group_names, "perm_id"))

            for (guild in c(target_guilds, "netprimprod")) {
                output_file <- ifelse(guild == "netprimprod", "ecosystem_indices_", "model_outputs_")
                output_extract_func <- ifelse(guild == "netprimprod", extract_ecosystem_indices_output, StrathE2E2:::extract_annual_domain_output)
                shap_effects <- e2e_driver_variance_analysis(permutation_plan, sub_output_path, guild, output_file_base = output_file, output_extraction_func = output_extract_func)
                shapley_results[shapley_results$output == guild, ]$shapley_effect <- shap_effects
            }

            write.csv(shapley_results, file.path(sub_output_path, "shapley_effects.csv"), row.names = FALSE)
        }
    }
}

# production_outputs <- c("netprimprod")
# for (model_parm in model_parameterisations) {
#     for (decade in decades) {
#         for (ssp in ssps) {
#             sub_output_path <- file.path(output_path, glue("{model_parm}_parameterisation"), ssp, decade)
#             cat("Assessing Shapley Effects for", model_parm, ssp, decade)

#             shapley_results <- expand.grid(variable_group = variable_group_names, output = production_outputs, stringsAsFactors = FALSE)
#             shapley_results$shapley_effect <- NA
#             permutation_plan <- read.csv(file.path(sub_output_path, "permutation_plan.csv"), col.names = c(variable_group_names, "perm_id"))

#             for (guild in production_outputs) {
#                 shap_effects <- e2e_driver_variance_analysis(permutation_plan, sub_output_path, guild, output_file_base = "ecosystem_indices_", output_extraction_func = extract_ecosystem_indices_output)
#                 shapley_results[shapley_results$output == guild, ]$shapley_effect <- shap_effects
#             }

#             write.csv(shapley_results, file.path(sub_output_path, "shapley_effects.csv"), row.names = FALSE)
#         }
#     }
# }
