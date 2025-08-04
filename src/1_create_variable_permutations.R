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

#### Build model permutations for the within decade - ESM/SSP variable swaps ----

# Build all 4 base variants per decade
variant_options <- expand.grid(ESM = ESMs, SSP = SSPs, stringsAsFactors = FALSE)
ssp_only_options <- data.frame(ESM = rep(ESMs[1], 2), SSP = SSPs)

# Create hybrid combinations within one decade: 4^6*2 = 8,192 runs
single_decade_combinations <- expand.grid(
    append(replicate(6, seq_len(nrow(variant_options)), simplify = FALSE), list(seq_len(nrow(ssp_only_options)))),
    stringsAsFactors = FALSE
)

single_decade_combinations <- as.data.frame(
    t(apply(single_decade_combinations, 1, function(row) {
        map_ESM_SSP_vars <- sapply(as.numeric(row[1:6]), function(x) map_variant_index(variant_options, x))
        map_SSP_var <- map_variant_index(ssp_only_options, as.numeric(row[7]))

        return(c(map_ESM_SSP_vars, map_SSP_var))
    }))
)
names(single_decade_combinations) <- blocks

hybrid_sources <- data.frame(matrix(
    NA,
    nrow = nrow(single_decade_combinations) * length(decades),
    ncol = length(blocks) + 1
))
names(hybrid_sources) <- c("perm_id", names(single_decade_combinations))
perm_id <- 1

for (decade in decades) {
    for (i in seq_len(nrow(single_decade_combinations))) {
        row <- single_decade_combinations[i, ]
        row <- sapply(row, function(x) paste0(decade, "-", x))

        hybrid_sources[perm_id, ] <- c("perm_id" = perm_id, row)
        perm_id <- perm_id + 1
    }
}
write.csv(hybrid_sources, "./outputs/within_decade_ssp_esm_permutations.csv")

#### Build model permutations for the within SSP-ESM variant - decadal variable swaps ----
decade_options <- data.frame(decade = decades)
esm_ssp_options <- expand.grid(ESM = ESMs, SSP = SSPs, stringsAsFactors = FALSE)
esm_ssp_options <- paste0(esm_ssp_options$ESM, "-", esm_ssp_options$SSP)

# Create hybrid combinations within one ESM-SSP combo: 3^7 = 2,187 runs
single_esm_ssp_combinations <- expand.grid(
    replicate(7, seq_len(nrow(decade_options)), simplify = FALSE),
    stringsAsFactors = FALSE
)

single_esm_ssp_combinations <- as.data.frame(
    t(apply(single_esm_ssp_combinations, 1, function(row) {
        sapply(as.numeric(row), function(x) map_variant_index(decade_options, x))
    }))
)
names(single_esm_ssp_combinations) <- blocks

hybrid_sources <- data.frame(matrix(
    NA,
    nrow = nrow(single_esm_ssp_combinations) * length(esm_ssp_options),
    ncol = length(blocks) + 1
))
names(hybrid_sources) <- c("perm_id", names(single_esm_ssp_combinations))
perm_id <- 1

for (esm_ssp in esm_ssp_options) {
    for (i in seq_len(nrow(single_esm_ssp_combinations))) {
        row <- single_esm_ssp_combinations[i, ]
        row <- sapply(row, function(x) paste0(x, "-", esm_ssp))

        hybrid_sources[perm_id, ] <- c("perm_id" = perm_id, row)
        perm_id <- perm_id + 1
    }
}
write.csv(hybrid_sources, "./outputs/within_esm_ssp_decade_permutations.csv")
