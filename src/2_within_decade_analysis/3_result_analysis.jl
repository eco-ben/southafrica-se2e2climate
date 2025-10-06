# using GlobalSensitivity
using CSV, DataFrames
# using Parquet
# using DataFrames
# using Combinatorics
using Statistics
using MultivariateStats
using MLDataUtils
using AlgebraOfGraphics
using GLMakie

within_decade_perms = CSV.read("../outputs/2010-2019_permutation_plan.csv", DataFrame)
rename!(
    within_decade_perms,
    "0" => "constant",
    "1" => "light",
    "2" => "temperature",
    "3" => "river_outputs",
    "4" => "vertical_mixing",
    "5" => "water_flows",
    "6" => "nutrient_concentrations",
    "7" => "atmospheric_nutrient_flux"
)
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
    "Cetaceans"
]

decade = "2010-2019"
shapley_effects = CSV.read("../outputs/2010-2019_shapley_effects.csv", DataFrame)

# Plot Shapley Effect values for different model outputs
output_plot = data(shapley_effects[shapley_effects.output .∈ [guilds[1:5]], :]) * mapping(:variable_group, :shapley_effect, row=:output) * visual(BarPlot, direction=:x)
fig = draw(output_plot)
fig

output_plot = data(shapley_effects[shapley_effects.output .∈ [guilds[6:13]], :]) * mapping(:variable_group, :shapley_effect, row=:output) * visual(BarPlot, direction=:x)
fig = draw(output_plot)
fig
