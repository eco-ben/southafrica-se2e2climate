library(ShapleyValue)
library(MASS)
library(tidyverse)
library(arrow)
source("./src/build_model_forcing_functions.R")

# Package example:
data <- Boston # Dataframe or matrix
y <- data$medv # Vector of outputs
x <- as.data.frame(data[, 5:8]) # dataframe or matrix of input values with rows for permutations and columns for variables
shapleyvalue(y, x)
shapleyvalue(y, x)[2, ] # Normalised shapley values

master_forcings <- load_all_possible_drivers("../../StrathE2E_workspace/Models/South_Africa_MA/", "South_Africa_MA")

find_encoded_val <- function(full_code_value, encoding_scheme, disregarded_pattern) {
    split_code_value <- gsub(disregarded_pattern, full_code_value, replacement = "")

    return(unname(encoding_scheme[split_code_value]))
}

summarise_annual_results <- function(df, variable, summarisation_method) {
    result <- summarisation_method(df[, variable])
    return(result)
}

# Analyse within decade permutations for each climate scenario level:
within_decade_perms <- read.csv("./outputs/within_decade_ssp_esm_permutations.csv")
input_variables <- colnames(within_decade_perms)[colnames(within_decade_perms) != "perm_id" & colnames(within_decade_perms) != "X"]
ssp_esm_encoding <- c("CNRM-ssp126" = 1, "GFDL-ssp126" = 2, "CNRM-ssp370" = 3, "GFDL-ssp370" = 4)

# For 2010-2019:
decade_perms_2010_2019 <- within_decade_perms[str_detect(within_decade_perms$light, "2010-2019"), ]
for (var in input_variables) {
    decade_perms_2010_2019[, var] <- sapply(
        decade_perms_2010_2019[, var],
        function(x) find_encoded_val(x, ssp_esm_encoding, "(\\d{4})-(\\d{4})-")
    )
}

test <- read_parquet("./outputs/within_decade_permutations/model_outputs_perm_1.parq")
annual_surface_phyt <- summarise_annual_results(test[test$Description == "Surface_layer_phytoplankton", ], "Model_annual_mean", function(x) x)

result_files <- lapply(1:200, function(i) result_files[[i]] <- str_glue("./outputs/within_decade_permutations/model_outputs_perm_{i}.parq"))
annual_surface_phyt <- sapply(result_files, function(x) {
    df <- read_parquet(x)
    summarise_annual_results(df[df$Description == "Surface_layer_phytoplankton", ], "Model_annual_mean", function(x) x)
})
annual_surface_phyt <- unlist(annual_surface_phyt)

shapleyvalue(annual_surface_phyt, decade_perms_2010_2019[1:200, ])





# Analyse within decade permutations for each climate scenario level:
within_esm_ssp_perms <- read.csv("./outputs/within_esm_ssp_decade_permutations.csv")
input_variables <- colnames(within_esm_ssp_perms)[colnames(within_esm_ssp_perms) != "perm_id" & colnames(within_esm_ssp_perms) != "X"]
decade_encoding <- c("2010-2019" = 1, "2030-2039" = 2, "2060-2069" = 3)

# For CNRM-SSP370:
esm_ssp_perms_CNRM_ssp370 <- within_esm_ssp_perms[str_detect(within_esm_ssp_perms$light, "CNRM-ssp126"), ]
for (var in input_variables) {
    esm_ssp_perms_CNRM_ssp370[, var] <- sapply(
        esm_ssp_perms_CNRM_ssp370[, var],
        function(x) find_encoded_val(x, decade_encoding, "-[[:upper:]]{4}-[[:lower:]]{3}\\d{3}")
    )
}

result_files <- list(1:200)
result_files <- lapply(1:200, function(i) result_files[[i]] <- str_glue("./outputs/across_decade_permutations/model_outputs_perm_{i}.parq"))
annual_surface_phyt <- sapply(result_files, function(x) {
    df <- read_parquet(x)
    summarise_annual_results(df[df$Description == "Surface_layer_phytoplankton", ], "Model_annual_mean", function(x) x)
})
annual_surface_phyt <- unlist(annual_surface_phyt)

shapleyvalue(annual_surface_phyt, esm_ssp_perms_CNRM_ssp370[1:200, ])

ggplot() +
    geom_point(data = master_forcings[master_forcings$variable %in% variable_groups["nutrient_conc"][[1]], ], aes(x = decade, y = value, color = SSP, size = ESM)) +
    facet_wrap(~variable)
