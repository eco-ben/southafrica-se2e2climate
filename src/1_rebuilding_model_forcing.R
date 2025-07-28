library(StrathE2E2)
library(tidyverse)
library(docstring)

# se2e_workspace <- "../../StrathE2E_workspace"
# analysis_workspace <- "C:/Users/kbb25108/OneDrive - University of Strathclyde/Documents/repos/southafrica-se2e2climate/"

# setwd(se2e_workspace)

se2e_physics_names <- c(
    "month", "sslight", "so_logespm", "si_logespm",
    "so_temp", "d_temp", "si_temp", "rivervol", "logkvert",
    "mixlscale", "d_so_upwelling", "so_d_downwelling", "so_inflow",
    "d_inflow", "si_inflow", "si_outflow", "so_si_flow",
    "s1_pdist", "s2_pdist", "s3_pdist", "d1_pdist", "d2_pdist",
    "d3_pdist", "Inshore_waveheight", "DO_logkvert", "DO_mixlscale",
    "DO_d_upwelling", "d_DO_downwelling"
)
se2e_chemistry_names <- c(
    "month",
    "so_nitrate", "so_ammonia", "so_phyt", "so_detritus",
    "d_nitrate", "d_ammonia", "d_phyt", "d_detritus",
    "si_nitrate", "si_ammonia", "si_phyt", "si_detritus",
    "rivnitrate", "rivammonia", "rivdetritus",
    "so_atmnitrate", "so_atmammonia",
    "si_atmnitrate", "si_atmammonia",
    "si_othernitrate", "si_otherammonia",
    "so_othernitrate", "so_otherammonia",
    "DO_nitrate", "DO_ammonia", "DO_detritus"
)

load_all_possible_drivers <- function(model_base_path, model_name, physics_names = se2e_physics_names, chemistry_names = se2e_chemistry_names) {
    #' Load all possible physics and chemistry driving values from model variants.
    #'
    #' @description This function creates a master dataframe containing all variants of
    #' physics and chemistry drivers. (Note that some drivers do not vary across model variants).
    #'
    #' @param model_base_path character. Base model folder that contains variant subfolders.
    #' @param model_name character. Model name contained in driving data file names.
    #' @param physics_names character. Correct names to apply to loaded physics file (copied from StrathE2E2).
    #' @param chemistry_names character. Correct names to apply to loaded chemistry file (copied from StrathE2E2).
    #' @return Longform dataframe containing all possible variant driving data for chemistry and physics variables.
    #'
    #' @references StrathE2E2
    #' @examples
    #' load_all_possible_drivers("Models/South_Africa_MA/", "South_Africa_MA")

    possible_variants <- dir(model_base_path)
    # decades <- c("2010-2019", "2030-2039", "2060-2069")
    # possible_variants <- possible_variants[
    #     str_detect(possible_variants, decades[1]) |
    #         str_detect(possible_variants, decades[2]) |
    #         str_detect(possible_variants, decades[3])
    # ]
    possible_variants <- possible_variants[str_detect(possible_variants, "-[:upper:]{4}-[:lower:]{3}\\d{3}")]
    decades <- unique(str_extract(possible_variants, "\\d{4}-\\d{4}"))
    ESMs <- unique(str_extract(possible_variants, "[:upper:]{4}"))
    SSPs <- unique(str_extract(possible_variants, "[:lower:]{3}\\d{3}"))

    master_forcings <- expand.grid(
        month = 1:12,
        decade = decades,
        SSP = SSPs,
        ESM = ESMs
    )

    for (variant in possible_variants) {
        ESM <- str_extract(variant, "[:upper:]{4}")
        SSP <- str_extract(variant, "[:lower:]{3}\\d{3}")
        decade <- str_extract(variant, "\\d{4}-\\d{4}")

        physics <- read.csv(paste0(
            model_base_path,
            decade, "-", ESM, "-", SSP,
            "/Driving/",
            "physics_", toupper(model_name), "_",
            decade, "-", ESM, "-", SSP,
            ".csv"
        ))
        names(physics) <- physics_names
        physics$decade <- decade
        physics$SSP <- SSP
        physics$ESM <- ESM

        chemistry <- read.csv(paste0(
            model_base_path,
            decade, "-", ESM, "-", SSP,
            "/Driving/",
            "chemistry_", toupper(model_name), "_",
            decade, "-", ESM, "-", SSP,
            ".csv"
        ))
        names(chemistry) <- chemistry_names
        chemistry$decade <- decade
        chemistry$SSP <- SSP
        chemistry$ESM <- ESM

        if (which(variant == possible_variants) == 1) {
            master_forcings <- left_join(master_forcings, physics, by = c("month", "decade", "SSP", "ESM"))
            master_forcings <- left_join(master_forcings, chemistry, by = c("month", "decade", "SSP", "ESM"))
        } else {
            master_forcings <- rows_update(master_forcings, physics, by = c("month", "decade", "SSP", "ESM"))
            master_forcings <- rows_update(master_forcings, chemistry, by = c("month", "decade", "SSP", "ESM"))
        }
    }

    master_forcings <- pivot_longer(
        master_forcings,
        cols = names(master_forcings)[!names(master_forcings) %in% c("month", "SSP", "ESM", "decade")],
        names_to = "variable",
        values_to = "value"
    )
}

rebuild_model_drivers <- function(template_model, master_forcing_values, variable_sources, physics_names = se2e_physics_names, chemistry_names = se2e_chemistry_names) {
    #' Rebuild model drivers from different sources.
    #'
    #' @description Rebuild template_model drivers using different variant values from variable_sources
    #'
    #' @param template_model StrathE2E2 model. Template e2e_read() model to recreate model structure.
    #' @param master_forcing_values dataframe. Master dataframe containing all possible variants of driving data.
    #' (col-names = month, decade, SSP, ESM, variable, value).
    #' @param variable_sources dataframe. Dataframe containing desired value sources for each driving variable.
    #' (col-names = variable, SSP, ESM, decade).
    #' @param physics_names character. Physics variable names (copied from StrathE2E2).
    #' @param chemistry_names character. Chemistry variable names (copied from StrathE2E2).
    #' @return StrathE2E2 model object containing driving data from variant sources indicated in variable_sources.
    #'
    #' @references StrathE2E2
    #' @examples
    #' rebuild_model_drivers(
    #'     e2e_read(model.name = "South_Africa_MA", model.variant = "2010-2019-CNRM-ssp126"),
    #'     load_all_possible_drivers("./Models/South_Africa_MA/", "South_Africa_MA"),
    #'     variable_sources
    #' )

    data <- StrathE2E2:::elt(template_model, "data")
    current_physics <- StrathE2E2:::elt(data, "physics.drivers")
    current_chemistry <- StrathE2E2:::elt(data, "chemistry.drivers")

    new_physics <- current_physics
    for (p_name in physics_names) {
        if (p_name == "month") {
            next
        }
        var_source <- variable_sources[variable_sources$variable == p_name, ]

        if (paste0(var_source$decade, "-", var_source$ESM, "-", var_source$SSP) == template_model$setup$model.variant) {
            next
        }

        new_vals <- master_forcings[
            master_forcings$decade == var_source$decade &
                master_forcings$ESM == var_source$ESM &
                master_forcings$SSP == var_source$SSP &
                master_forcings$variable == p_name,
        ]

        new_physics[, p_name] <- new_vals[sort(new_vals$month), ]$value
    }

    new_chemistry <- current_chemistry
    for (c_name in chemistry_names) {
        if (c_name == "month") {
            next
        }
        var_source <- variable_sources[variable_sources$variable == c_name, ]

        if (paste0(var_source$decade, "-", var_source$ESM, "-", var_source$SSP) == template_model$setup$model.variant) {
            next
        }

        new_vals <- master_forcings[
            master_forcings$decade == var_source$decade &
                master_forcings$ESM == var_source$ESM &
                master_forcings$SSP == var_source$SSP &
                master_forcings$variable == c_name,
        ]
        new_physics[, c_name] <- new_vals[sort(new_vals$month), ]$value
    }

    new_data <- data
    new_data[["physics.drivers"]] <- new_physics
    new_data[["chemistry.drivers"]] <- new_chemistry

    new_model <- list(
        setup = template_model$setup,
        data = new_data
    )

    return(new_model)
}

# master <- load_all_possible_drivers("./Models/South_Africa_MA/", "South_Africa_MA")
# unique_variables <- unique(master$variable)
# variable_sources <- data.frame(
#     variable = unique_variables,
#     SSP = rep("ssp126", length(unique_variables)),
#     ESM = rep("CNRM", length(unique_variables)),
#     decade = rep("2060-2069", length(unique_variables))
# )

# t_model <- e2e_read(model.name = "South_Africa_MA", model.variant = possible_variants[1], models.path = "Models")
# new_model <- rebuild_model_drivers(t_model, master, variable_sources)
