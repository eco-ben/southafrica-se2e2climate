library(StrathE2E2)
library(tidyverse)
library(furrr)

# n_workers <- 20
# plan(multisession, workers = n_workers)

source("./project_config.R")

variants <- list.dirs(
    str_glue("{source("./project_config.R")}S_Benguela_MA/"),
    recursive = FALSE
)
variants <- variants[!str_detect(variants, "2010-2015")] # Remove 2010-2015 variants
variants <- variants[str_detect(variants, "ssp")] # Include only variant folders with ESM/SSP info (remove 2010-2019 base variant)

variants <- str_split_i(variants, str_glue("{source("./project_config.R")}S_Benguela_MA/"), 2)
results_path_base <- "./outputs/initial_runs/"
sapply(variants, function(x) dir.create(file.path(results_path_base, x)))

npp_mc_post <- function(
    result_dir,
    ident_baseline,
    ident_change,
    esm_ssp,
    creds = c(0.025, 0.25, 0.5, 0.75, 0.975)
) {
    baseline_whole <- list.files(
        file.path(result_dir, ident_baseline, "CredInt"),
        pattern = ident_baseline,
        full.names = TRUE
    )
    change_whole <- list.files(
        file.path(result_dir, ident_change, "CredInt"),
        pattern = ident_change,
        full.names = TRUE
    )

    baseline_whole <- read.csv(baseline_whole[grep(
        baseline_whole,
        pattern = "CredInt_cumulative_wholeannualflux"
    )])
    change_whole <- read.csv(change_whole[grep(
        change_whole,
        pattern = "CredInt_cumulative_wholeannualflux"
    )])

    prop_change_whole <- baseline_whole[, c("iteration", "likelihood")]
    prop_change_whole$netprimprod <- (change_whole[,
        "Phytoplankton_net_primary_production"
    ] -
        baseline_whole[, "Phytoplankton_net_primary_production"]) /
        baseline_whole[, "Phytoplankton_net_primary_production"]

    results <- data.frame(rep(0, length(creds) + 1))
    rownames(results) <- c(
        "maxlik",
        "lowlimit",
        "lowquart",
        "median",
        "uppquart",
        "upplimit"
    )

    results[, 1] <- rep(0, length(creds) + 1)
    colnames(results)[1] <- "netprimprod"
    results[1, 1] <- prop_change_whole$netprimprod[1]

    results$netprimprod[2:nrow(results)] <- StrathE2E2:::GetCredInt(
        prop_change_whole$netprimprod,
        prop_change_whole$likelihood,
        creds,
        var = "netprimprod",
        plotgraph = FALSE
    )

    results <- cbind(data.frame("cred_int" = rownames(results)), results)
    results$esm_ssp <- esm_ssp
    results$ident_baseline <- ident_baseline
    results$ident_changed <- ident_change

    write.csv(
        results,
        file.path(
            result_dir,
            "mass_files",
            paste0(
                "CredInt_processed_whole_npp_change_from_",
                ident_baseline,
                "_to_",
                ident_change,
                ".csv"
            )
        ),
        row.names = FALSE
    )
}

for (esm_ssp in esm_ssps) {
    scen_variants <- variants[str_detect(variants, esm_ssp)]
    baseline <- scen_variants[str_detect(scen_variants, "2010-2019")]
    scen_variants <- scen_variants[scen_variants != baseline]

    walk(scen_variants, function(scen_var) {
        npp_mc_post(
            result_dir = "./outputs/initial_runs/monte_carlo_results/S_Benguela_MA/",
            baseline,
            scen_var,
            esm_ssp
        )
    })

    walk(variants[str_detect(variants, esm_ssp)], function(variant) {
        var_files <- list.files(
            file.path(result_dir, variant, "CredInt"),
            pattern = variant,
            full.names = TRUE
        )
        flux <- read.csv(var_files[grep(
            var_files,
            pattern = "CredInt_cumulative_wholeannualflux"
        )])
        results <- data.frame(rep(0, length(creds) + 1))
        rownames(results) <- c(
            "maxlik",
            "lowlimit",
            "lowquart",
            "median",
            "uppquart",
            "upplimit"
        )

        results[, 1] <- rep(0, length(creds) + 1)
        colnames(results)[1] <- "netprimprod"
        results[1, 1] <- flux[, "Phytoplankton_net_primary_production"][1]

        results$netprimprod[2:nrow(results)] <- StrathE2E2:::GetCredInt(
            flux[, "Phytoplankton_net_primary_production"],
            flux$likelihood,
            c(0.025, 0.25, 0.5, 0.75, 0.975),
            var = "netprimprod",
            plotgraph = FALSE
        )
        results <- cbind(data.frame("cred_int" = rownames(results)), results)
        write.csv(
            results,
            file.path(
                result_dir,
                "mass_files",
                paste0(
                    "CredInt_processed_whole_npp_",
                    variant,
                    ".csv"
                )
            ),
            row.names = FALSE
        )
    })
}
