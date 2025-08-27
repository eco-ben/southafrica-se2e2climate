source("src/model_permutations_common.R")

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
