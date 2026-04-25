# southafrica-se2e2climate

[![DOI](https://zenodo.org/badge/1021318483.svg)](https://doi.org/10.5281/zenodo.19751528)

Analysis code for StrathE2E2 climate change paper using the Southern Benguela
Mission Atlantic implementation.

# Setup

Full recreation of this project requires that users have R v4.5.2, Julia v1.11.6
and Quarto installed. Set up the Julia project in the usual Julian way by 
instantiating the project.

Start Julia from the project root (the main project directory, where the README file sits).

```shell
julia --project=.
```

Followed by:

```julia
] instantiate
```

The above will setup the project and install all dependencies.

Note that the specific version of Rasters.jl and Makie.jl currently has an error preventing
precompilation of the Makie extension for Rasters.jl (RastersMakieExt.jl). An error message
will appear, but can be safely ignored as we do not rely on this extension for any
visualizations.

All scripts are expected to be run from the `src` directory.
The current directory can be changed without leaving Julia by hitting semi-colon.

```julia
; cd src
```

## Structure

Repository structure:
``` code
southafrica-se2e2climate/
├─ src/     # Analysis scripts
├─ outputs/     # StrathE2E model and analysis results 
├─ data/    # Input data (that is not produced by code in `src`) but is required for plotting
├─ paper/   # Reproducible quarto manuscript via the rendering of .qmd file
├─ figs/   # Output directory for figures*
├─ project_config.R    # Configuration file containing the path to the model file directory
├─ .gitignore
├─ Project.toml     # Julia project spec
├─ Manifest.toml    # Julia project spec
└─ README.md    # this file
```

\* Three figures have been produced externally and have been uploaded to the repository
in the `/figs/` directory. "ssp126_daily_masked_smoothed_twilight_shifted_r_difference_multi.png",
"ssp370_daily_masked_smoothed_twilight_shifted_r_difference_multi.png" and 
"StrathE2E_paper_method_outline.png".

`project_config.R` File must be modified to contain the path (relative from this project root)
to the parent folder containing the Southern Benguela StrathE2E implementation files
(e.g. "../../Model_files/", where Model_files contains the folder "S_Benguela_MA/").

## Required data files

The StrathE2E Southern Benguela model implementation files are required to be stored in
a separate directory. These files can be downloaded from 

The following required data files are not created as StrathE2E outputs and are
needed in the `/data/` folder for plotting.

- `GEBCO_2020.nc` : GEBCO 2020 grid bathymetry data available from https://www.gebco.net/data-products/gridded-bathymetry-data/gebco-2020
- `fleet-daily.parq` : Global Fishing Watch apparent fishing effort data collated for the South Africa fleets into daily point data. File available from https://doi.org/10.5281/zenodo.19750821
- `Habitats.gpkg` : StrathE2E habitat GeoPackage file created from sediment maps, available from https://doi.org/10.5281/zenodo.19750821
- `ne_10m_land shapefile` : Natural Earth 10m resolution land shapefile available from https://www.naturalearthdata.com/downloads/10m-physical-vectors/
- `sau_landings_strath_gears.parq` : Sea Around Us data for the region collated into StrathE2E gears and guilds. File available from https://doi.org/10.5281/zenodo.19750821
- `distance-from-shore.tif` : GeoTiff file containing values for the distance from ocean cells to the South African land border. File available from https://doi.org/10.5281/zenodo.19750821

# Running analyses

R scripts are primarily used to run analyses using the StrathE2E2 package and associated
functions, this includes Monte Carlo analyses, sensitivity analysis of demersal fish
and Shapley Effect calculations.

Julia scripts are primarily used to run analyses on model output data and create
figures for publication. Julia is used for it's analysis speed and easy ability
to create and customise professional figures.

The scripts in the project should be run in the following order:
`1_initial_climate_change/`
- `1_initial_climate_change_montecarlo.R` : Monte Carlo simulations for all variants of the model.
- `1_2_monte_carlo_npp.R` : Custom collation of Net Primary Production data after Monte Carlo simulations
- `1_2_demfish_sens.R` : Demersal fish sensitivity analyses

`2_across_ESM_analysis/`
- `1_permutation_runs.R` : Generating and running the climate change sensitivity permutations across NEMO-ERSEM forcing models.

`3_across_decade_analysis/`
- `1_permutation_runs.R` : Generating and running the climate change sensitivity permutations across decades.

`1_initial_climate_change/2_initial_cc_assessment.jl` : Analyses and figures for initial climate change results, including selection of guilds for later analyses
`2_across_ESM_analysis/3_result_analysis.jl` : Shapley Effect figures for the guilds identified as important across NEMO-ERSEM forcing models.
`3_across_decade_analysis/3_result_analysis.jl` : Figures created using permutations across decades.
`4_publication_plots.jl` : Additional plotting required for publication (creation of supplementary parameterisation plots).
`5_plotting_environmental_drivers.jl` : Supplementary plots for environmental driving data.

Additional script `analysis_common.jl` contains common elements used across all Julia scripts for analysis and figure creation.

**Important note:** For running scripts `1_initial_climate_change_montecarlo.R`
and `1_2_monte_carlo_npp.R` the StrathE2E2 package should use the GitLab repository
branch `phyt_mass_mcprocess`. For running scripts named `1_permutation_runs.R` the
StrathE2E2 package should use the GitLab repository branch `climate-analysis`.