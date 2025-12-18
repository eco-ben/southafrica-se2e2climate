# southafrica-se2e2climate
Analysis code for StrathE2E2 climate change paper using the Southern Benguela - South Africa MA implementation.

# Setup

Full recreation of this project requires that users have Julia v1.11.6 and Quarto installed.
Set up the project in the usual Julian way by instantiating the project.

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
ADRIA-reef-indicators/
├─ src/     # Analysis scripts
├─ outputs/     # ADRIA results sets and analysis / figure outputs
├─ data/
├─ paper/   # Reproducible quarto manuscript via the rendering of .qmd file
├─ .gitignore
├─ config.toml  # configuration file for project
├─ Project.toml     # Julia project spec
├─ Manifest.toml    # Julia project spec
└─ README.md    # this file
```

The ADRIA Domain data package `GBR_2024_10_15_HighResCoralStress` is required to run ADRIA
for this project.

The path to this domain folder should be specified in the config.toml file.

The following structure and files should be maintained in the outputs, figs and data folders:
``` code
ADRIA-reef-indicators
│───outputs
│   └───ADRIA_results
│       └───HighResCoralStress
│           ├───processed_model_outputs     # Folder containing scenario-median timeseries arrays for each GCM
│           │       ├───median_cover_ACCESS-CM2.nc
│           │       ├───median_cover_ACCESS-ESM1-5.nc
│           │       ├───median_cover_EC-Earth3-Veg.nc
│           │       ├───median_cover_GFDL-CM4.nc
│           │       └───median_cover_NorESM2-MM.nc
│           ├───HighResCoralStress_ADRIA_scens_ACCESS-CM2.csv       # CSV files containing the ADRIA-CoralBlox input parameter scenarios used in 1_.jl for each GCM
│           ├───HighResCoralStress_ADRIA_scens_ACCESS-ESM1-5.csv
│           ├───HighResCoralStress_ADRIA_scens_EC-Earth3-Veg.csv
│           ├───HighResCoralStress_ADRIA_scens_GFDL-CM4.csv
│           ├───HighResCoralStress_ADRIA_scens_NorESM2-MM.csv
│           ├───clustered_reefs_carbonate.gpkg      # Reef spatial data with cluster timeseries information attached (before context layers are attached in 3_.jl).
│           └───analysis_context_layers_carbonate.gpkg      # Geopackage containing the required timeseries clustering, connectivity and carbonate budget values for each reef
│───figs    # Directory directly containing general GCM-wide and paper-methods figures, as well as GCM subdirectories
│   ├───ACCESS-CM2
│   │   ├───bioregion
│   │   ├───gbr
│   │   └───management_area
│   ├───ACCESS-ESM1-5
│   │   ├───bioregion
│   │   ├───gbr
│   │   └───management_area
│   ├───EC-Earth3-Veg
│   │   ├───bioregion
│   │   ├───gbr
│   │   └───management_area
│   ├───GFDL-CM4
│   │   ├───bioregion
│   │   ├───gbr
│   │   └───management_area
│   └───NorESM2-MM
│       ├───bioregion
│       ├───gbr
│       └───management_area
└───data
    ├─ GBRMPA_Management_Areas.gpkg    # File containing polygons of GBRMPA management areas
    ├─ GBRMPA_Reef_Features.gpkg    # File containing polygons of GBRMPA reef feature data
    └─ GBRMPA_Reefal_Bioregions.gpkg   # File containing polygons of GBRMPA reefal bioregions.
```
Additionally, if available, ADRIA Result Sets should be located in the outputs/ADRIA_Results/HighResCoralStress/ folder.
