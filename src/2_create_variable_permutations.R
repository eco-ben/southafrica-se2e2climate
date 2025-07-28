source("src/rebuilding_model_forcing.R")

master <- load_all_possible_drivers("../../StrathE2E_workspace/Models/South_Africa_MA/", "South_Africa_MA")

decades <- c("2010-2019", "2030-2039", "2060-2069") # Define decades that I want to investigate (currently only 2010-2019, 2030-39 and 2060-69)
ESMs <- unique(master$ESM)
SSPs <- unique(master$SSP)
blocks <- c(
    "light", "temperature", "river_outputs", "vertical_mixing", "water_flows", "nutrient_conc", "other_nut_flux", "atm_nut_flux"
)

# Build all 4 base variants per decade
variant_options <- expand.grid(ESM = ESMs, SSP = SSPs, stringsAsFactors = FALSE)
ssp_only_options <- data.frame(ESM = rep(ESMs[1], 2), SSP = SSPs)

# Create hybrid combinations within one decade: 4^7 = 16,384
single_decade_combinations <- expand.grid(
    append(replicate(7, seq_len(nrow(variant_options)), simplify = FALSE), list(seq_len(length(SSPs)))),
    stringsAsFactors = FALSE
)

# Function to map numeric index to actual ESM/SSP pair
map_variant_index <- function(options, index) {
    paste(options[index, ], collapse = "-")
}

single_decade_combinations <- as.data.frame(
    t(apply(single_decade_combinations, 1, function(row) {
        map_ESM_SSP_vars <- sapply(as.numeric(row[1:7]), function(x) map_variant_index(variant_options, x))
        map_SSP_var <- map_variant_index(ssp_only_options, as.numeric(row[8]))

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

write.csv(hybrid_sources, "./outputs/climate_variant_permutations.csv")
