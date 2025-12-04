using Colors
using CairoMakie
using Statistics
import Random

Random.seed!(1234)

# Define constant guilds and guild individual colours
guilds = [
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
    "Cetaceans",
    "netprimprod",
    "Demersal_fish_larvae"
]
guild_clean_names = Dict(zip(guilds, replace.(guilds, "_" => " ")))
guild_clean_names["netprimprod"] = "Net Primary Production"
guild_individual_colours = Dict(zip(guilds, distinguishable_colors(length(guilds))))

flux_guilds = Dict(
    "Surface_layer_phytoplankton" => "phyt",
    "Deep_layer_phytoplankton" => "phyt",
    "Omnivorous_zooplankton" => "omnivzoo",
    "Carnivorous_zooplankton" => "carnzoo",
    "Benthos_susp/dep_feeders" => "benths",
    "Benthos_carn/scav_feeders" => "benthc",
    "Planktivorous_fish" => "pfish",
    "Migratory_fish" => "mfish",
    "Demersal_fish" => "dfish",
    "Birds" => "bird",
    "Pinnipeds" => "seal",
    "Cetaceans" => "ceta",
    "Demersal_fish_larvae" => "dfishlar"
)

# Define constant variables for importance analysis
variable_levels = [
    "0" => "constant",
    "1" => "light",
    "2" => "temperature",
    "3" => "river_outputs",
    "4" => "vertical_mixing",
    "5" => "boundary_flows",
    "6" => "nutrient_concentrations",
    "7" => "atmospheric_nutrient_flux"
]
variables = ["light", "river_outputs", "nutrient_concentrations", "atmospheric_nutrient_flux", "temperature", "vertical_mixing", "boundary_flows"]
variable_clean_names = Dict(zip(variables, replace.(variables, "_" => " ")))
variable_colours = Dict(zip(variables, distinguishable_colors(length(variables))))

variable_groups = Dict(
    "light" => ["sslight"], 
    "river_outputs" => ["rivervol", "rivnitrate", "rivammonia"], 
    "nutrient_concentrations" => ["so_nitrate", "so_ammonia", "so_phyt", "so_detritus", "d_nitrate", "d_ammonia", "d_phyt", "d_detritus", "si_nitrate", "si_ammonia", "si_phyt", "si_detritus"],
    "atmospheric_nutrient_flux" => ["so_atmnitrate", "so_atmammonia", "si_atmnitrate", "si_atmammonia"], 
    "temperature" => ["so_temp", "d_temp", "si_temp"], 
    "vertical_mixing" => ["logkvert", "d_so_upwelling", "so_d_downwelling"], 
    "boundary_flows" => ["so_inflow", "d_inflow", "si_inflow", "si_outflow", "so_si_flow"]
)

# Define constant ESM/SSP colours
ESM_categories = ["GFDL" => "GFDL-ESM4", "CNRM" => "CNRM-CM6-1-HR"]
ESM_colors = ["GFDL" => Makie.wong_colors()[1], "CNRM" => Makie.wong_colors()[2]]
SSP_categories = ["ssp126" => "SSP1-2.6", "ssp370" => "SSP3-7.0"]
SSP_linestyles = ["ssp126" => :solid, "ssp370" => :dash]

ESMs = ["CNRM", "GFDL"]
SSPs = ["ssp126", "ssp370"]

ESM_SSPs = ["CNRM-ssp126", "CNRM-ssp370", "GFDL-ssp126", "GFDL-ssp370"]
ESM_SSP_categories = [
    "CNRM-ssp126" => "CNRM-CM6-1-HR SSP1-2.6", 
    "CNRM-ssp370" => "CNRM-CM6-1-HR SSP3-7.0", 
    "GFDL-ssp126" => "GFDL-ESM4 SSP1-2.6", 
    "GFDL-ssp370" => "GFDL-ESM4 SSP3-7.0"
]
ESM_SSP_colors = [
    "CNRM-ssp126" => Makie.wong_colors()[1], 
    "CNRM-ssp370" => Makie.wong_colors()[2], 
    "GFDL-ssp126" => Makie.wong_colors()[3], 
    "GFDL-ssp370" => Makie.wong_colors()[4]
]

fontsize = 7
dpi = 300

# Size of 1cm in pixels relative to 1 CSS px, see:
# - https://docs.makie.org/dev/how-to/match-figure-size-font-sizes-and-dpi
# - https://docs.makie.org/dev/explanations/figure#Figure-size-and-resolution
centimetre = 37.7952755906

# Convert fontsize to pixel measurement
pt = 1.33 # size of 1 pt in pixels
fontsize = fontsize * pt

inch = 96 # size of 1 inch in pixels
dpi = dpi / inch

function rescale!(x)
    if eltype(x) <: Integer
        x = convert.(Float64, x)
    end

    μ = mean(x)
    σ = std(x)

    # avoid divide-by-zero — MLDataUtils behavior
    σ == 0 && (σ = 1.0)

    x .= (x .- μ) ./ σ

    return μ, σ
end