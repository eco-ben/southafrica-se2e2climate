# using GlobalSensitivity
using CSV, DataFrames
# using Parquet
# using DataFrames
# using Combinatorics
using Statistics
using MultivariateStats
using AlgebraOfGraphics
using GLMakie

include("../analysis_common.jl")

output_path = "../outputs/across_esm_ssp_permutations"
figs_path = "../figs/across_esm_ssp_permutations"
decades = ["2010-2019", "2030-2039", "2060-2069"]

for decade in decades
    
    sub_output_path = joinpath(output_path, decade)
    sub_figs_path = joinpath(figs_path, decade)

    permutation_plan = CSV.read(joinpath(sub_output_path, "permutation_plan.csv"), DataFrame)
    rename!(
        permutation_plan,
        variable_levels
    )

    files = readdir(sub_output_path; join=true)
    files = files[contains.(files, ["model_outputs_"])]
    perm_outputs = [CSV.read(files[contains.(files, ["perm_$(pid).csv"])], DataFrame) for pid in permutation_plan.perm_id]

    perm_plan_outputs = permutation_plan
    for guild in guilds
        output = [out_df[out_df.Description .== [guild], :Model_annual_mean] for out_df in perm_outputs]
        output = vcat(output...)
        perm_plan_outputs[!, guild] = output
    end

    shapley_effects = CSV.read(joinpath(sub_output_path, "shapley_effects.csv"), DataFrame)

    shapley_effects_wide = unstack(shapley_effects, :output, :variable_group, :shapley_effect)
    rename!(shapley_effects_wide, "water_flows" => "boundary_flows")

    dist_mat = pairwise(Cityblock(), Matrix(shapley_effects_wide[:, names(shapley_effects_wide) .∉ [["output", "cluster"]]])')
    # dist_mat = dist_mat./ maximum(dist_mat)

    # shapley_effects_wide.cluster = kmeans(Matrix{Float64}(Matrix(shapley_effects_wide[:, Not(:output, :cluster)])'), 5).assignments

    shapley_effects_wide.cluster = kmedoids(dist_mat, 4).assignments

    # function kluster_auto_k(X, clust_func; kmin=2, kmax=10, qindex=:silhouettes)
    #     k_range = kmin:kmax
    #     n_k = length(k_range)
    #     res_x, k_x, score_x = [], [], [] 

    #     for (i, k) in enumerate(k_range)
    #         res = clust_func(X, k)
            
    #         s = mean(clustering_quality(X, res.assignments; quality_index=qindex))
    #         # if s > best_score
    #         #     best_k, best_score, best_res = k, s, res
    #         # end
    #         push!(res_x, res)
    #         push!(k_x, k)
    #         push!(score_x, s)
    #     end

    #     return res_x, k_x, score_x
    # end

    # shapley_effects_wide.cluster = first(kluster_auto_k(dist_mat, kmedoids; qindex=:silhouettes)).assignments
    shapley_effects_wide = sort(shapley_effects_wide, [:cluster, :output])
    shapley_plot = shapley_effect_parallel_plot(shapley_effects_wide)
    save(joinpath(sub_figs_path, "guild_clustered_importances.png"), shapley_plot, px_per_unit = dpi)

end

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
shapley_effects = CSV.read("../outputs/across_esm_ssp_permutations/2010-2019/shapley_effects.csv", DataFrame)

# Plot Shapley Effect values for different model outputs
output_plot = data(shapley_effects[shapley_effects.output .∈ [guilds], :]) * mapping(:variable_group, :shapley_effect, row=:output) * visual(BarPlot, direction=:x)
fig = draw(output_plot)
fig

output_plot = data(shapley_effects[shapley_effects.output .∈ [guilds[6:13]], :]) * mapping(:variable_group, :shapley_effect, row=:output) * visual(BarPlot, direction=:x)
fig = draw(output_plot)
fig
