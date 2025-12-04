"""
This script is for the analysis and plotting of the across ESM permutation results. These
analyses are focussed on the guilds identified as being important for across ESM ecosystem
state variation in the initial PCA analysis (2_initial_cc_assessment.jl).
"""

using CSV, DataFrames
using Statistics
using MultivariateStats
using AlgebraOfGraphics
using CairoMakie
using Clustering
using CategoricalArrays
using Distances

include("../analysis_common.jl")

output_path = "../outputs/across_esm_permutations"
figs_path = "../figs/across_esm_permutations"

all_decades = ["2010-2019", "2020-2029", "2030-2039", "2040-2049", "2050-2059", "2060-2069"]

function variant_shapley_effects(esm, ssp, dec)
    shapley_effects = joinpath(output_path, "$(esm)_parameterisation", ssp, dec, "shapley_effects.csv")
    shapley_effects = CSV.read(shapley_effects, DataFrame)
    rename!(shapley_effects, ["output" => "guild", "variable_group" => "variable"])
    shapley_effects.variable = ifelse.(shapley_effects.variable .== "water_flows", "boundary_flows", shapley_effects.variable)

    shapley_effects.ESM_param .= esm
    shapley_effects.SSP .= ssp
    shapley_effects.decade .= dec

    return shapley_effects
end

all_shapley_effects = vcat([
    variant_shapley_effects(esm, ssp, dec) for (esm, ssp, dec) in Iterators.product(ESMs, SSPs, all_decades)
]...)
all_shapley_effects.SSP_decade = all_shapley_effects.SSP .* "-" .* all_shapley_effects.decade
all_shapley_effects = all_shapley_effects[all_shapley_effects.variable .!= "constant", :]

aggregated_shap = all_shapley_effects[all_shapley_effects.guild .∈ [esm_sep_guilds; "netprimprod"], :]
aggregated_shap = combine(
    groupby(aggregated_shap, [:guild, :ESM_param, :variable]), 
    :shapley_effect => mean,
    :shapley_effect => minimum,
    :shapley_effect => maximum,
    :shapley_effect => (x -> quantile(x, 0.75)) => :shapley_effect_75,
    :shapley_effect => (x -> quantile(x, 0.25)) => :shapley_effect_25
)
aggregated_shap = sort(aggregated_shap, :shapley_effect_mean)
aggregated_shap.variable_clean_name = getindex.([variable_clean_names], aggregated_shap.variable)
aggregated_shap.guild_clean_name = getindex.([guild_clean_names], aggregated_shap.guild)
aggregated_shap.jitter = rand(-0.2:0.05:0.2, nrow(aggregated_shap))

all_shapley_effects 

unique_variables = unique(aggregated_shap.variable_clean_name)
aggregated_shap.x_position = Vector{Int64}(indexin(aggregated_shap.variable_clean_name, unique_variables))
aggregated_shap.std_bar_line = [
    [[r.x_position + r.jitter, r.x_position + r.jitter],
    [max(0, r.shapley_effect_mean - r.shapley_effect_minimum), min(1, r.shapley_effect_mean + r.shapley_effect_maximum)]]
for r in eachrow(aggregated_shap)]

esm_guilds = unique(aggregated_shap.guild)
esm_guild_colours = Dict(zip(esm_guilds, Makie.wong_colors()[eachindex(esm_guilds)]))

cnrm = aggregated_shap[aggregated_shap.ESM_param .== "CNRM", :]
gfdl = aggregated_shap[aggregated_shap.ESM_param .== "GFDL", :]

fig = Figure(fontsize = fontsize, size = (18.42centimetre, 12centimetre))

ax1 = Axis(fig[1,1], title = "CNRM parameterisation", ylabel = "mean Shapley effect", xticks = (eachindex(unique_variables), unique_variables), xticklabelrotation=π/4)
scatter!(ax1, cnrm.x_position .+ cnrm.jitter, cnrm.shapley_effect_mean, color = getindex.([esm_guild_colours], cnrm.guild))
map(x -> lines!(ax1, cnrm.std_bar_line[x][1], cnrm.std_bar_line[x][2], color = getindex(esm_guild_colours, cnrm.guild[x])), eachindex(eachrow(cnrm)))

ax2 = Axis(fig[1,2], title = "GFDL parameterisation", xticks = (eachindex(unique_variables), unique_variables), xticklabelrotation=π/4)
scatter!(ax2, gfdl.x_position .+ gfdl.jitter, gfdl.shapley_effect_mean, color = getindex.([esm_guild_colours], gfdl.guild))
map(x -> lines!(ax2, gfdl.std_bar_line[x][1], gfdl.std_bar_line[x][2], color = getindex(esm_guild_colours, gfdl.guild[x])), eachindex(eachrow(gfdl)))

esm_guild_elements = [[
    MarkerElement(color=esm_guild_colours[guild], marker=:circle), 
    LineElement(color=esm_guild_colours[guild])
] for guild in esm_guilds]
Legend(fig[2, :], esm_guild_elements, getindex.([guild_clean_names], esm_guilds), nbanks=2, orientation=:horizontal)

save("../figs/across_esm_permutations/mean_shapley_effects.png", fig, px_per_unit=dpi)
