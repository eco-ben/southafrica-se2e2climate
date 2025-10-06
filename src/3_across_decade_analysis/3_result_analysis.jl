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

esm_ssp = "CNRM-ssp126"
shapley_effects = CSV.read("../outputs/CNRM-ssp126_shapley_effects.csv", DataFrame)

# Plot Shapley Effect values for different model outputs
output_plot = data(shapley_effects[shapley_effects.output .∈ [guilds[1:5]], :]) * mapping(:variable_group, :shapley_effect, row=:output) * visual(BarPlot, direction=:x)
fig = draw(output_plot)
fig

output_plot = data(shapley_effects[shapley_effects.output .∈ [guilds[6:13]], :]) * mapping(:variable_group, :shapley_effect, row=:output) * visual(BarPlot, direction=:x)
fig = draw(output_plot)
fig

# Perform PCA on Shapley Effect values with outputs as samples to identify clusters of outputs based on variable importance
shapley_effects_wide = unstack(shapley_effects, :output, :variable_group, :shapley_effect)

M_pca = fit(PCA, Matrix{Float64}(Matrix(shapley_effects_wide[:, Not(:output)])'); maxoutdim=2)
y_pca = predict(M_pca, Matrix{Float64}(Matrix(shapley_effects_wide[:, Not(:output)])'))
pca_val_df = DataFrame(output = shapley_effects_wide.output, PC1 = y_pca[1, :], PC2 = y_pca[2, :])

pca_plot = data(pca_val_df) * mapping(:PC1, :PC2, color=:output) * visual(Scatter)
fig = draw(pca_plot) 

## Test applying kmeans clustering to the guild outputs based on variable importance values (clustering on PCA values also results in the same clusters)
test = kmeans(y_pca, 3)
test = kmeans(Matrix{Float64}(Matrix(shapley_effects_wide[:, Not(:output)])'), 3)
pca_val_df.cluster = categorical(test.assignments)

pca_plot = data(pca_val_df) * mapping(:PC1, :PC2, color=:cluster) * visual(Scatter)
fig = draw(pca_plot)

legend = only(filter(x -> x isa Legend, fig.content))
