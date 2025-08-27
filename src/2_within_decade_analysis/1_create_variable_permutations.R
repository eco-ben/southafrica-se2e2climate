source("src/model_permutations_common.R")

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
