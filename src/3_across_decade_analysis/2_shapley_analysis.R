library(StrathE2E2)
library(tidyverse)
library(glue)

models_path <- "../../StrathE2E_workspace/Models/"
output_path <- "./outputs/across_decade_permutations"

model <- e2e_read("South_Africa_MA", "2010-2019-CNRM-ssp126", models.path = models_path)
ESM_SSP <- c("CNRM-ssp126", "CNRM-ssp370", "GFDL-ssp126", "GFDL-ssp370")
variable_group_names <- c(
    "constant",
    "light",
    "temperature",
    "river_outputs",
    "vertical_mixing",
    "water_flows",
    "nutrient_concentrations",
    "atmospheric_nutrient_flux"
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

# Add Shapley effects for net primary production
extract_ecosystem_indices_output <- function(csv_file, output_id) {
  result <- read.csv(csv_file)
  result <- result[result$Description == output_id, "Value"]
  
  return(result)
}

for (esm_ssp in ESM_SSP) {
    message("Assessing Shapley effects for ", esm_ssp)

    sub_output_path <- paste0(output_path, "/", esm_ssp, "/")
    shapley_results <- expand.grid(variable_group = variable_group_names, output = c(target_guilds, "netprimprod"), stringsAsFactors = FALSE)
    shapley_results$shapley_effect <- NA
    permutation_plan <- read.csv(paste0(sub_output_path, "permutation_plan.csv"), col.names = c(variable_group_names, "perm_id"))

    for (guild in c(target_guilds, "netprimprod")) {
      output_file <- ifelse(guild == "netprimprod", "ecosystem_indices_", "model_outputs_")
      output_extract_func <- ifelse(guild == "netprimprod", extract_ecosystem_indices_output, StrathE2E2:::extract_annual_domain_output)
      shap_effects <- e2e_driver_variance_analysis(permutation_plan, sub_output_path, guild, output_file_base = output_file, output_extraction_func = output_extract_func)
      shapley_results[shapley_results$output == guild, ]$shapley_effect <- shap_effects
    }

    write.csv(shapley_results, paste0(sub_output_path, "shapley_effects.csv"), row.names = FALSE)
}
