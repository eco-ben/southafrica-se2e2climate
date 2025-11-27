# using GlobalSensitivity
using CSV, DataFrames
# using Parquet
# using DataFrames
# using Combinatorics
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

scat = data(all_shapley_effects[all_shapley_effects.guild .∈ [pca_loadings[pca_loadings.similarity .>= cosd(30), :guild]], :]) * mapping(:variable, :shapley_effect, color=:SSP_decade, marker=:guild, layout=:ESM_param) * visual(Scatter)
fig = draw(scat)

test = all_shapley_effects[all_shapley_effects.guild .∈ [pca_loadings[pca_loadings.esm_sep_similarity .>= cosd(30), :guild]], :]
test = combine(
    groupby(test, [:guild, :ESM_param, :variable]), 
    :shapley_effect => median,
    :shapley_effect => minimum,
    :shapley_effect => maximum,
    :shapley_effect => (x -> quantile(x, 0.75)) => :shapley_effect_75,
    :shapley_effect => (x -> quantile(x, 0.25)) => :shapley_effect_25
)
test = sort(test, :shapley_effect_median)
test.variable_clean_name = getindex.([variable_clean_names], test.variable)
test.guild_clean_name = getindex.([guild_clean_names], test.guild)
test.jitter = rand(-0.2:0.05:0.2, nrow(test))

unique_variables = unique(test.variable_clean_name)
test.x_position = Vector{Int64}(indexin(test.variable_clean_name, unique_variables))
test.std_bar_line = [
    [[r.x_position + r.jitter, r.x_position + r.jitter],
    [max(0, r.shapley_effect_25), min(1, r.shapley_effect_75)]]
for r in eachrow(test)]

esm_guilds = unique(test.guild)
esm_guild_colours = Dict(zip(esm_guilds, Makie.wong_colors()[eachindex(esm_guilds)]))

cnrm = test[test.ESM_param .== "CNRM", :]
gfdl = test[test.ESM_param .== "GFDL", :]

fig = Figure(fontsize = fontsize, size = (18.42centimetre, 12centimetre))

ax1 = Axis(
    fig[1,1], 
    title = "CNRM parameterisation", 
    ylabel = "median Shapley effect", 
    xticks = (eachindex(unique_variables), unique_variables), 
    xticklabelrotation=π/4,
    ylabelpadding=15
)
scatter!(ax1, cnrm.x_position .+ cnrm.jitter, cnrm.shapley_effect_median, color = getindex.([esm_guild_colours], cnrm.guild))
map(x -> lines!(ax1, cnrm.std_bar_line[x][1], cnrm.std_bar_line[x][2], color = getindex(esm_guild_colours, cnrm.guild[x])), eachindex(eachrow(cnrm)))

ax2 = Axis(
    fig[1,2], 
    title = 
    "GFDL parameterisation", 
    xticks = (eachindex(unique_variables), 
    unique_variables), xticklabelrotation=π/4
)
scatter!(ax2, gfdl.x_position .+ gfdl.jitter, gfdl.shapley_effect_median, color = getindex.([esm_guild_colours], gfdl.guild))
map(x -> lines!(ax2, gfdl.std_bar_line[x][1], gfdl.std_bar_line[x][2], color = getindex(esm_guild_colours, gfdl.guild[x])), eachindex(eachrow(gfdl)))

esm_guild_elements = [[
    MarkerElement(color=esm_guild_colours[guild], marker=:circle), 
    LineElement(color=esm_guild_colours[guild])
] for guild in esm_guilds]
Legend(fig[2, :], esm_guild_elements, getindex.([guild_clean_names], esm_guilds), nbanks=2, orientation=:horizontal)
save("../figs/across_esm_permutations/mean_shapley_effects.png", fig, px_per_unit=dpi)


fig_opts = (;
    fontsize=fontsize, size = (18.42centimetre, 13centimetre)
)
ax_opts = (;
    xticklabelrotation=π/4
)
legend_opts = (;
    position=:bottom,
    nbanks=2
)
scat = data(test) * mapping(:guild_clean_name, :shapley_effect_median, color=:variable, layout=:ESM_param) * visual(Scatter)
fig = draw(scat; figure=fig_opts, axis=ax_opts, legend=legend_opts)



test = all_shapley_effects[all_shapley_effects.guild .∈ [pca_loadings[pca_loadings.similarity .>= cosd(30), :guild]], :]
test = combine(
    groupby(test, [:guild, :ESM_param, :variable]), 
    :shapley_effect => mean,
    :shapley_effect => minimum,
    :shapley_effect => maximum,
    :shapley_effect => (x -> quantile(x, 0.75)) => :shapley_effect_75,
    :shapley_effect => (x -> quantile(x, 0.25)) => :shapley_effect_25
)
test = sort(test, :shapley_effect_mean)
test.variable_clean_name = getindex.([variable_clean_names], test.variable)
test.guild_clean_name = getindex.([guild_clean_names], test.guild)
test.jitter = rand(-0.2:0.05:0.2, nrow(test))

all_shapley_effects 

unique_variables = unique(test.variable_clean_name)
test.x_position = Vector{Int64}(indexin(test.variable_clean_name, unique_variables))
test.std_bar_line = [
    [[r.x_position + r.jitter, r.x_position + r.jitter],
    [max(0, r.shapley_effect_mean - r.shapley_effect_minimum), min(1, r.shapley_effect_mean + r.shapley_effect_maximum)]]
for r in eachrow(test)]

esm_guilds = unique(test.guild)
esm_guild_colours = Dict(zip(esm_guilds, Makie.wong_colors()[eachindex(esm_guilds)]))

cnrm = test[test.ESM_param .== "CNRM", :]
gfdl = test[test.ESM_param .== "GFDL", :]

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

# fig_opts = (;
#     fontsize = fontsize,
#     size = (18.42centimetre, 14centimetre)
# )
# scale = scales(
#     Y = (; label = "Variable groups"),
#     X = (; label = "Shapley Effect"),
#     Color = (; label = "Guilds")
# )
# axis_opts = (; aspect=1)
# legend_opts = (; position=:bottom, tellheight=false, tellwidth=false, nbanks=2)
# scat = data(test[test.guild .∈ [pca_loadings[pca_loadings.similarity .>= cosd(30), :guild]], :]) * mapping(:shapley_effect_mean, :variable_clean_name, color=:guild_clean_name, layout=:ESM_param) * visual(Scatter)
# fig = draw(scat; legend=(;position=:bottom, nbanks=2))

for esm_ssp in ESM_SSPs
    
    sub_output_path = joinpath(output_path, esm_ssp)
    sub_figs_path = joinpath(figs_path, esm_ssp)

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

    
    base_variant = "2010-2019-" * esm_ssp
    comp_variant = "2060-2069-" * esm_ssp

    var_permutation_outputs = vcat([
        collate_variable_interactions(permutation_plan, perm_plan_outputs, base_variant, comp_variant, var, variables)
        for var in variables
    ]...)
    
    for guild in guilds
        figure = plot_guild_var_interactions(var_permutation_outputs, variables, guild)
        guild = replace(guild, "/" => "-")
        save(joinpath(sub_figs_path, "$(guild)_variable_interactions.png"), figure, px_per_unit=dpi)
    end
end


function shapley_effect_parallel_plot(
    shapley_effects_wide
)
    
    clusters = shapley_effects_wide.cluster 
    cluster_colours = Dict(zip(unique(clusters), Makie.wong_colors()[2:length(unique(clusters))+1]))

    fig = Figure(size = (18.42centimetre, 13.4centimetre), fontsize=fontsize)
    # variables = ["light", "temperature", "river_outputs", "vertical_mixing", "water_flows", "nutrient_concentrations", "atmospheric_nutrient_flux"]
    
    ax1 = Axis(fig[1,1], xticks=(eachindex(variables), variables), ylabel = "Shapley Effects", xticklabelrotation=π/8)
    # ax2 = Axis(fig[2,1], xticks=(eachindex(variables), variables), ylabel = "Standardised Shapley Effects")

    guild_lines = [[(rand(-0.05:0.01:0.05) + v, r[var]) for (v, var) in enumerate(variables)] for r in eachrow(shapley_effects_wide)]

    for (v, variable) in enumerate(variables)
        jitter = first.(getindex.(guild_lines, v))
        # scatter!(ax1, jitter, shapley_effects_wide[:, variable], color=getindex.([guild_colours], shapley_effects_wide.output), markersize=20)
        scatter!(ax1, jitter, shapley_effects_wide[:, variable], color=:grey, markersize=5)
        # scatter!(ax2, jitter, shapley_effects_standard[:, variable], color=getindex.([guild_colours], shapley_effects_standard.output), markersize=10, alpha=0.7)
    end
    map(x -> lines!(ax1, guild_lines[x], color=getindex(cluster_colours, shapley_effects_wide.cluster[x])), eachindex(clusters))

    # guild_elements = [MarkerElement(color=guild_colours[guild], marker=:circle) for guild in shapley_effects_wide.output]
    line_elements = [PolyElement(color=cluster_colours[cluster]) for cluster in shapley_effects_wide.cluster]
    # Legend(fig[1, 2], guild_elements, shapley_effects_wide.output)
    Legend(fig[1, 1], line_elements, getindex.([guild_clean_names], shapley_effects_wide.output), "Guild Clusters", halign=:left, valign=:top, tellwidth=false, tellheight=false)

    return fig
end

# Plot Shapley Effect values for different model outputs
output_plot = data(shapley_effects[shapley_effects.output .∈ [guilds], :]) * mapping(:variable_group, :shapley_effect, layout=:output) * visual(BarPlot, direction=:x)
fig = draw(output_plot)
fig

save("../figs/example_shapley_plots.png", fig)

constant_perms = findall([all(row) for row in eachrow(permutation_plan.light .== permutation_plan[:, Not(:constant, :perm_id)])])
constant_perm_results = [CSV.read(joinpath(sub_output_path, "model_outputs_perm_$(i).csv"), DataFrame) for i in constant_perms]
for (i, df) in enumerate(constant_perm_results)
    perm_id = constant_perms[i]
    perm_desc = permutation_plan[perm_id, :light]
    df.perm_id .= constant_perms[i]
    df.perm_desc .= perm_desc
end
constant_perm_results = vcat(constant_perm_results...)

function plot_guild_shap_change(
    shapley_effects_df, 
    constant_perm_res_df, 
    title; 
    biomass_output_col = "Model_annual_mean", 
    shap_col = "shapley_effect", 
    biomass_ylabel = "Biomass [mMNpergWW]", 
    shap_ylabel = "Shapley Effect",
    biomass_xcol = "perm_desc",
    shap_xcol = "variable_group"
)
    fig = Figure(size = (1000, 400))

    ax1 = Axis(fig[1,1], xticks = (1:nrow(constant_perm_res_df), constant_perm_res_df[:, biomass_xcol]), ylabel = biomass_ylabel)
    barplot!(ax1, 1:nrow(constant_perm_res_df), constant_perm_res_df[:, biomass_output_col])

    ax2 = Axis(fig[2,1], xticks = (1:nrow(shapley_effects_df), shapley_effects_df[:, shap_xcol]), ylabel = shap_ylabel, xticklabelrotation = 2/π)
    barplot!(ax2, 1:nrow(shapley_effects_df), shapley_effects_df[:, shap_col])

    Label(fig[0, :], title, tellwidth = false, tellheight = false)
    rowsize!(fig.layout, 0, Relative(0.01))

    return fig
end

function plot_guild_shap_var(
    shapley_effects_df, 
    constant_perm_res_df,
    perm_plan_outputs,
    guild, 
    title; 
    biomass_output_col = "Model_annual_mean", 
    shap_col = "shapley_effect", 
    biomass_ylabel = "Biomass [mMNpergWW]", 
    shap_ylabel = "Shapley Effect",
    biomass_xcol = "perm_desc",
    shap_xcol = "variable_group"
)
    variance_lines = [
        [(1, minimum(perm_plan_outputs[:, guild])), (1, maximum(perm_plan_outputs[:, guild]))],
        [(0.9999, minimum(perm_plan_outputs[:, guild])), (1.0001, minimum(perm_plan_outputs[:, guild]))],
        [(0.9999, maximum(perm_plan_outputs[:, guild])), (1.0001, maximum(perm_plan_outputs[:, guild]))]
    ]
    colours = Dict(zip(constant_perm_res_df[:, biomass_xcol], Makie.wong_colors()[1:nrow(constant_perm_res_df)]))
    elements = [MarkerElement(color = colours[var], marker=:circle) for var in constant_perm_res_df[:, biomass_xcol]]

    fig = Figure(size = (25centimetre, 15centimetre), fontsize = fontsize)
    ax1 = Axis(fig[1,1]; ylabel = biomass_ylabel, xticksvisible=false, xticklabelsvisible=false)
    lines!.([ax1], variance_lines; color=:black, linewidth=3)
    scatter!(ax1, fill(1, nrow(constant_perm_res_df)), constant_perm_res_df[:, biomass_output_col], color = getindex.([colours], constant_perm_res_df[:, biomass_xcol]), markersize=20)

    Legend(fig[2, 1], elements, constant_perm_res_df[:, biomass_xcol], orientation = :horizontal, nbanks = nrow(constant_perm_res_df))

    ax2 = Axis(fig[1,2], xticks = (1:nrow(shapley_effects_df), shapley_effects_df[:, shap_xcol]), ylabel = shap_ylabel)
    barplot!(ax2, 1:nrow(shapley_effects_df), shapley_effects_df[:, shap_col])

    Label(fig[0, :], title, tellwidth = false, tellheight = false)
    rowsize!(fig.layout, 0, Relative(0.01))
    colsize!(fig.layout, 1, Relative(0.1))
    rowsize!(fig.layout, 2, Relative(0.1))

    return fig
end

for guild in guilds
    clean_name = replace(guild, "_" => " ")
    figure = plot_guild_shap_var(
        shapley_effects[shapley_effects.output .== guild, :],
        constant_perm_results[constant_perm_results.Description .== guild, :],
        perm_plan_outputs,
        guild,
        clean_name
    )
    guild = replace(guild, "/" => "-")
    save("../figs/across_decade_permutations/CNRM-ssp370/$(guild)_shapley.png", figure)
end

output_plot = data(shapley_effects[shapley_effects.output .∈ [guilds[6:13]], :]) * mapping(:variable_group, :shapley_effect, row=:output) * visual(BarPlot, direction=:x)
fig = draw(output_plot)
fig

# Perform PCA on Shapley Effect values with outputs as samples to identify clusters of outputs based on variable importance
shapley_effects_wide = unstack(shapley_effects, :output, :variable_group, :shapley_effect)

shapley_effects_standard = shapley_effects_wide
rescaled_vars = names(shapley_effects_wide[:, Not(:output, :cluster)])
for col in rescaled_vars
    if all(shapley_effects_standard[:, col] .== 0.0) continue end
    μ, σ = rescale!(shapley_effects_standard[!, col])
end

M_pca = fit(PCA, Matrix{Float64}(Matrix(shapley_effects_wide[:, Not(:output)])'); maxoutdim=8)
y_pca = predict(M_pca, Matrix{Float64}(Matrix(shapley_effects_wide[:, Not(:output)])'))
pca_val_df = DataFrame(output = shapley_effects_wide.output, PC1 = y_pca[1, :], PC2 = y_pca[2, :], PC3 = y_pca[3, :])

colours = Dict(zip(pca_val_df.output, distinguishable_colors(length(pca_val_df.output))))
elements = [MarkerElement(color=colours[guild], marker=:circle) for guild in pca_val_df.output]

fig = Figure(size = (800, 800), fontsize = fontsize)
ax = Axis(fig[1,1])
scatter!(ax, pca_val_df.PC1, pca_val_df.PC2, color=getindex.([colours], pca_val_df.output))
Legend(fig[1,2], elements, pca_val_df.output, nbanks=1)

test_fig = data(shapley_effects) * mapping(:variable_group, :shapley_effect, color=:output) * visual(Scatter)
test_fig = draw(test_fig)

save("../figs/across_decade_permutations/CNRM-ssp370/guild_shapley_pca.png", fig)

## Test applying kmeans clustering to the guild outputs based on variable importance values (clustering on PCA values also results in the same clusters)
test = kmeans(y_pca, 3)
test = kmeans(Matrix{Float64}(Matrix(shapley_effects_wide[:, Not(:output)])'), 4)
pca_val_df.cluster = categorical(test.assignments)

pca_plot = data(pca_val_df) * mapping(:PC1, :PC2, :PC3, text=:output => verbatim, color=:cluster) * visual(Makie.Text; fontsize = 10) 
fig = draw(pca_plot, figure = (; fontsize = fontsize))

colours = Dict(zip(unique(pca_val_df.cluster), distinguishable_colors(length(unique(pca_val_df.cluster)))))
fig = Figure()
ax1 = Axis(fig[1,1], xlabel="PC1", ylabel="PC2")
text!(ax1, pca_val_df.PC1, pca_val_df.PC2, text=pca_val_df.output, color=getindex.([colours], pca_val_df.cluster))
ax2 = Axis(fig[1,2], xlabel="PC1", ylabel="PC3")
text!(ax2, pca_val_df.PC1, pca_val_df.PC3, text=pca_val_df.output, color=getindex.([colours], pca_val_df.cluster))
ax3 = Axis(fig[2,1], xlabel="PC2", ylabel="PC3")
text!(ax3, pca_val_df.PC2, pca_val_df.PC3, text=pca_val_df.output, color=getindex.([colours], pca_val_df.cluster))

# legend_entries = [MarkerElement(color = Makie.wong_colors()[3:8][cl], marker=:circle) for cl in pca_val_df.cluster.refs]
# Legend(fig.figure[1,2], legend_entries, pca_val_df.output)
# Legend(fig.figure[2, 1], [MarkerElement(color=col, marker=:circle) for col in Makie.wong_colors()[3:5]], ["Cluster 1", "Cluster 2", "Cluster 3"], orientation = :horizontal)
# rowsize!(fig.figure.layout, 2, Relative(0.05))

L = loadings(M_pca)
# ev = principalvars(M_pca) ./ sum(principalvars(M_pca))  # variance explained
# corr_circle = L .* sqrt.(ev[1:2]')                      # scale loadings
corr_circle = L
# scaling_factor = maximum(abs, y_pca)  # heuristic for visibility
# L_scaled = corr_circle .* scaling_factor
# corr_circle = L_scaled

fig = fig.figure
# ax2 = Axis(fig[2, 1], xlabel="PC1", ylabel="PC2")

# draw unit circle (correlation circle)
# θ = range(0, 2π, length=200)
# lines!(ax, cos.(θ), sin.(θ), color=:gray, linestyle=:dash)
annual_outputs = unique(shapley_effects.variable_group)
annual_outputs_colors = Dict(zip(annual_outputs, distinguishable_colors(length(annual_outputs))))

ax = Axis(
    fig.layout[2,1],
    limits = (extrema(y_pca[1, :]), extrema(y_pca[2, :]))
)
# Draw arrows for variables
for (i, output) in enumerate(annual_outputs)
    arrows2d!(ax, [mean(y_pca[1, :])], [mean(y_pca[2, :])], [corr_circle[i,1]], [corr_circle[i,2]], 
            shaftwidth=2, color=annual_outputs_colors[output])
    text!(ax, corr_circle[i,1], corr_circle[i,2], text=annual_outputs[i], align=(:left, :bottom), alpha=0.5)
end

Label(fig.layout[1,1, TopLeft()], "A", font=:bold)
Label(fig.layout[2,1, TopLeft()], "B", font=:bold)

save("../figs/across_decade_permutations/CNRM-ssp370/guild_clustered_shapley_pca.png", fig)

baseline_decade = "CNRM-ssp126"
decade_2 = "CNRM-ssp370"
variable = "temperature"
other_variables = names(within_decade_perms[:, Not([:perm_id, Symbol(variable)])])

baseline_decades = contains.(within_decade_perms[:, Not(:perm_id)], baseline_decade)
perm_all_baseline = within_decade_perms[[all(contains.(r, baseline_decade)) for r in eachrow(within_decade_perms)], :]
perm_var_changed = within_decade_perms[[(all(r))], :]

guilds = shapley_effects_wide.output
guild_colours = Dict(zip(unique(guilds), distinguishable_colors(length(unique(guilds)))))
clusters = pca_val_df.cluster
cluster_colours = Dict(zip(unique(clusters), distinguishable_colors(length(unique(clusters)))))

dat = shapley_effects_wide[:, Not([:output, :constant])]
n_vars = size(dat, 2)
angles = range(0, 2π, length = n_vars + 1)
radar_data = [vcat(collect(dat[i, :]), dat[i, 1]) for i in 1:size(dat, 1)]
var_labels = names(dat)

fig = Figure()
ax1 = PolarAxis(fig[1,1], thetaticks = (collect(angles)[Not(end)], var_labels), title="Coloured by guild")
ax2 = PolarAxis(fig[2,1], thetaticks = (collect(angles)[Not(end)], var_labels), title="Coloured by cluster")

for i in 1:size(dat, 1)
    lines!(ax1, angles, radar_data[i], color = guild_colours[guilds[i]], linewidth=2, alpha=0.7)
    
    lines!(ax2, angles, radar_data[i], color = cluster_colours[clusters[i]], linewidth=2, alpha=0.7)
    # Optional: add labels at the last point
    # text!(ax, angles[end], radar_data[i][end], labels[i], align = (:left, :center), fontsize=10)
end

Legend(fig[1,2], [PolyElement(color=col) for col in getindex.([guild_colours], guilds)], guilds)
Legend(fig[2,2], [PolyElement(color=col) for col in getindex.([cluster_colours], clusters)], guilds)

test_data = shapley_effects_wide
test_data.cluster = pca_val_df.cluster
test_data = sort(test_data, :cluster)

test_data_mat = test_data[:, Not([:output, :constant, :cluster])]

fig = Figure(size=(20centimetre,15centimetre), fontsize=fontsize+2)
ax1 = Axis(fig[1,1], yticks = (1:nrow(test_data_mat), test_data.output), xticks = (1:ncol(test_data_mat), names(test_data_mat)), limits=((0.5, n_vars+0.5), (0.5, nrow(test_data)+0.5)))
hm = heatmap!(ax1, Matrix(test_data_mat)')
ax2 = Axis(fig[1, 2], xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false, limits=(nothing, (0.5, nrow(test_data)+0.5)))
map(x -> hspan!(ax2, x-0.5, x+0.5, color = cluster_colours[test_data[x, :cluster]]), 1:nrow(test_data))
linkaxes!(ax1, ax2)
Colorbar(fig[1,3], hm, label="Standardised Shapley Effects")
colsize!(fig.layout, 2, Relative(0.1))
colgap!(fig.layout, 1, Relative(0))

cluster_colours = Dict(zip(unique(clusters), Makie.wong_colors()[2:length(unique(clusters))+1]))
test_data_pca = sort(pca_val_df, :cluster)
pca_lines = [[(1, r.PC1), (2, r.PC2), (3, r.PC3)] for r in eachrow(test_data_pca)]
guild_elements = [MarkerElement(color=guild_colours[guild], marker=:circle) for guild in test_data_pca.output]
cluster_elements = [MarkerElement(color=cluster_colours[cluster], marker=:circle) for cluster in unique(test_data_pca.cluster)]
#
fig = Figure(size = (30centimetre, 17centimetre), fontsize=fontsize)
ax = Axis(fig[1,1], xticks=(1:3, ["PC1", "PC2", "PC3"]))
# text!(ax, fill(1, nrow(test_data_pca)), test_data_pca.PC1, text = test_data_pca.output, align=(:center, :bottom))
# text!(ax, fill(2, nrow(test_data_pca)), test_data_pca.PC3, text = test_data_pca.output, align=(:center, :bottom))
# text!(ax, fill(3, nrow(test_data_pca)), test_data_pca.PC3, text = test_data_pca.output, align=(:center, :bottom))
map(x -> lines!(ax, pca_lines[x], color=getindex(cluster_colours, test_data_pca[x, :].cluster), alpha=0.6), eachindex(pca_lines))
scatter!(ax, fill(1, nrow(test_data_pca)), test_data_pca.PC1, color=getindex.([guild_colours], test_data_pca.output))
scatter!(ax, fill(2, nrow(test_data_pca)), test_data_pca.PC2, color=getindex.([guild_colours], test_data_pca.output))
scatter!(ax, fill(3, nrow(test_data_pca)), test_data_pca.PC3, color=getindex.([guild_colours], test_data_pca.output))

# gl = GridLayout(fig[1,2])
Legend(fig[1,2], guild_elements, test_data_pca.output)
Legend(fig[1,3], cluster_elements, ["Cluster $(i)" for i in unique(test_data_pca.cluster)])
colgap!(fig.layout, 2, Relative(0.005))

ax = Axis(
    fig.layout[2,1],
    limits = (extrema(y_pca[1, :]), extrema(y_pca[2, :]))
)
# Draw arrows for variables
for (i, output) in enumerate(annual_outputs)
    arrows2d!(ax, [mean(y_pca[1, :])], [mean(y_pca[2, :])], [corr_circle[i,1]], [corr_circle[i,2]], 
            shaftwidth=2, color=annual_outputs_colors[output])
    text!(ax, corr_circle[i,1], corr_circle[i,2], text=annual_outputs[i], align=(:left, :bottom), alpha=0.5)
end

save("../figs/across_decade_permutations/CNRM-ssp370/guild_clustered_importances.png", fig, px_per_unit=dpi)


# text!(ax, fill(1, nrow(test_data)), test_data.light, text = test_data.output, align=(:center, :bottom))
# text!(ax, fill(2, nrow(test_data)), test_data.temperature, text = test_data.output, align=(:center, :bottom))
# text!(ax, fill(3, nrow(test_data)), test_data.river_outputs, text = test_data.output, align=(:center, :bottom))
# text!(ax, fill(4, nrow(test_data)), test_data.vertical_mixing, text = test_data.output, align=(:center, :bottom))
# text!(ax, fill(5, nrow(test_data)), test_data.water_flows, text = test_data.output, align=(:center, :bottom))
# text!(ax, fill(6, nrow(test_data)), test_data.nutrient_concentrations, text = test_data.output, align=(:center, :bottom))
# text!(ax, fill(7, nrow(test_data)), test_data.atmospheric_nutrient_flux, text = test_data.output, align=(:center, :bottom))

# scatter!(ax, fill(1, nrow(test_data)), test_data.light, color=getindex.([guild_colours], test_data.output))
# scatter!(ax, fill(2, nrow(test_data)), test_data.temperature, color=getindex.([guild_colours], test_data.output))
# scatter!(ax, fill(3, nrow(test_data)), test_data.river_outputs, color=getindex.([guild_colours], test_data.output))
# scatter!(ax, fill(4, nrow(test_data)), test_data.vertical_mixing, color=getindex.([guild_colours], test_data.output))
# scatter!(ax, fill(5, nrow(test_data)), test_data.water_flows, color=getindex.([guild_colours], test_data.output))
# scatter!(ax, fill(6, nrow(test_data)), test_data.nutrient_concentrations, color=getindex.([guild_colours], test_data.output))
# scatter!(ax, fill(7, nrow(test_data)), test_data.atmospheric_nutrient_flux, color=getindex.([guild_colours], test_data.output))

# test = data(shapley_effects[shapley_effects.variable_group .!= "constant",:]) * mapping(:variable_group, :shapley_effect, color=:output) * visual(Scatter)
# draw(test)

function find_permutation(permutation_plan, base_variant, comparison_variant, var_vector, variables; permutation_id_col=:perm_id)
    if var_vector == "baseline"
        bool_vec = [all(collect(row) .== base_variant) for row in eachrow(permutation_plan[:, variables])]
    else
        other_variables = setdiff(variables, var_vector)
        bool_vec = [
            all(collect(row[var_vector]) .== comparison_variant) & 
            all(collect(row[other_variables]) .== base_variant) 
            for row in eachrow(permutation_plan[:, variables])
        ]
    end

    if all(.!bool_vec) 
        throw("Error finding permutation for variables $(var_vector)") 
    end

    return first(permutation_plan[bool_vec, permutation_id_col])
end

function collate_variable_interactions(permutation_plan, perm_plan_outputs, base_variant, comp_variant, variable, variables; permutation__id_col=:perm_id)
    others = setdiff(variables, [variable])
    interactions = vcat.(variable, others)
    
    variable_comparison_levels = DataFrame(variable_changed = vcat("baseline", [[variable]], interactions))
    variable_comparison_levels.perm_id = map(
        x -> find_permutation(permutation_plan, base_variant, comp_variant, x, variables),
        variable_comparison_levels.variable_changed
    )

    for guild in guilds
        variable_comparison_levels[!, guild] = perm_plan_outputs[variable_comparison_levels.perm_id, guild]
        variable_comparison_levels[!, "$(guild)_percent_change"] = 100 .- (variable_comparison_levels[variable_comparison_levels.variable_changed .== "baseline", guild] ./ 
            variable_comparison_levels[:, guild] .* 100)
    end

    variable_comparison_levels.labels = ifelse.(
        variable_comparison_levels.variable_changed .== "baseline",
        variable_comparison_levels.variable_changed,
        join.(variable_comparison_levels.variable_changed, ":")
    )
    variable_comparison_levels.analysis_variable .= variable

    return variable_comparison_levels
end

variables
variable = ["temperature"]



function plot_guild_var_interactions(permutation_outputs, variables, guild; guild_clean_names=guild_clean_names, title=guild_clean_names[guild])
    fig = Figure(size=(20centimetre, 20centimetre), fontsize=6pt)

    nrows = ceil(Int, sqrt(length(variables)))
    ncols = ceil(Int, length(variables) / nrows)
    indx = fldmod1.(1:nrows*ncols, ncols)

    for (v, variable) in enumerate(variables)
        variable_comparison_levels = permutation_outputs[permutation_outputs.analysis_variable .== variable, :]
        variable_comparison_levels = variable_comparison_levels[variable_comparison_levels.variable_changed .!= "baseline", :]

        ax = Axis(
            fig[indx[v]...], 
            xticks = (
                1:nrow(variable_comparison_levels), 
                replace.(variable_comparison_levels.labels, "_" => " ")
            ),
            title = var_clean_names[variable],
            xticklabelrotation = π / 7
        )
        barplot!(
            ax, 
            1:nrow(variable_comparison_levels), 
            variable_comparison_levels[:, "$(guild)_percent_change"],
            color=ifelse.(variable_comparison_levels[:, "$(guild)_percent_change"] .< 0, :red, :green),
            alpha=0.7
        )

        linkaxes!(filter(x -> x isa Axis, fig.content)...)
    end

    Label(fig[1:nrows, 0], "Percentage change from baseline\n permutation [%]", tellwidth=false, tellheight=false, rotation=pi/2)
    colsize!(fig.layout, 0, Relative(0.05))
    Label(fig[0, 1:ncols], title, font=:bold, tellwidth=false, tellheight=false)
    rowsize!(fig.layout, 0, Relative(0.05))

    return fig 
end

plot_guild_var_interactions(var_permutation_outputs, variables, "Planktivorous_fish")

ESM_SSPs

shapley_effects_all = [CSV.read(joinpath(output_path, esm_ssp, "shapley_effects.csv"), DataFrame) for esm_ssp in ESM_SSPs]
for (e, df) in enumerate(shapley_effects_all) df[!, "ESM_SSP"] .= ESM_SSPs[e] end
shapley_effects_all = vcat(shapley_effects_all...)

shapley_wide = unstack(shapley_effects_all[shapley_effects_all.variable_group .!= "constant", :], [:output, :ESM_SSP], :variable_group, :shapley_effect)

M_pca = fit(PCA, Matrix{Float64}(shapley_wide[:, Not(:output, :ESM_SSP)])', maxoutdim=2)
y_pca = predict(M_pca, Matrix{Float64}(shapley_wide[:, Not(:output, :ESM_SSP)])')
shapley_wide.PC1 = y_pca[1, :]
shapley_wide.PC2 = y_pca[2, :]

lines = data(shapley_wide) * mapping(:PC1, :PC2, color=:output) * (visual(Lines))
scatter = data(shapley_wide) * mapping(:PC1, :PC2, color=:output, marker=:ESM_SSP) * (visual(Scatter))

fig = draw(scatter; figure=fig_opts, legend=legend_opts)


function getellipsepoints(cx, cy, rx, ry, θ)
	t = range(0, 2*pi, length=100)
	ellipse_x_r = @. rx * cos(t)
	ellipse_y_r = @. ry * sin(t)
	R = [cos(θ) sin(θ); -sin(θ) cos(θ)]
	r_ellipse = [ellipse_x_r ellipse_y_r] * R
	x = @. cx + r_ellipse[:,1]
	y = @. cy + r_ellipse[:,2]
	(x,y)
end
function getellipsepoints(μ, Σ, confidence=0.95)
	quant = quantile(Chisq(2), confidence) |> sqrt
	cx = μ[1]
	cy =  μ[2]
	
	egvs = eigvals(Σ)
	if egvs[1] > egvs[2]
		idxmax = 1
		largestegv = egvs[1]
		smallesttegv = egvs[2]
	else
		idxmax = 2
		largestegv = egvs[2]
		smallesttegv = egvs[1]
	end

	rx = quant*sqrt(largestegv)
	ry = quant*sqrt(smallesttegv)
	
	eigvecmax = eigvecs(Σ)[:,idxmax]
	θ = atan(eigvecmax[2]/eigvecmax[1])
 	if θ < 0
		θ += 2*π
	end

	getellipsepoints(cx, cy, rx, ry, θ)
end

ESM_SSP_μ = [
    [
        mean(shapley_wide[shapley_wide.ESM_SSP .== esm_ssp, :PC1]), 
        mean(shapley_wide[shapley_wide.ESM_SSP .== esm_ssp, :PC2])
    ]
    for esm_ssp in ESM_SSPs
]
ESM_SSP_cov = [
    cov(Matrix(shapley_wide[shapley_wide.ESM_SSP .== esm_ssp, [:PC1, :PC2]]))
    for esm_ssp in ESM_SSPs
]
ESM_SSP_ellipse = getellipsepoints.(ESM_SSP_μ, ESM_SSP_cov, 0.999)
ESM_SSP_ellipse = [GB.Polygon(Point.(tuple.(points[1], points[2]))) for points in ESM_SSP_ellipse]
ESM_SSP_ellipse = DataFrame(ESM_SSP = ESM_SSPs, ellipse = ESM_SSP_ellipse)
ESM_SSP_convex = [
    GO.convex_hull(Point.(tuple.(
        shapley_wide[shapley_wide.ESM_SSP .== esm_ssp, :PC1], 
        shapley_wide[shapley_wide.ESM_SSP .== esm_ssp, :PC2]
    )))
    for esm_ssp in ESM_SSPs
]

dist_mat = pairwise(Cityblock(), Matrix(shapley_wide[:, variables])')
shapley_wide.cluster = kmedoids(dist_mat, 4).assignments
clusters = unique(shapley_wide.cluster)
cluster_colours = Dict(zip(clusters, Makie.wong_colors()[1:length(clusters)]))

cluster_μ = [
    [
        mean(shapley_wide[shapley_wide.cluster .== cluster, :PC1]), 
        mean(shapley_wide[shapley_wide.cluster .== cluster, :PC2])
    ]
    for cluster in clusters
]
cluster_cov = [
    cov(Matrix(shapley_wide[shapley_wide.cluster .== cluster, [:PC1, :PC2]]))
    for cluster in clusters
]
cluster_ellipse = getellipsepoints.(cluster_μ, cluster_cov, 0.999)
cluster_ellipse = [GB.Polygon(Point.(tuple.(points[1], points[2]))) for points in cluster_ellipse]
cluster_ellipse = DataFrame(cluster = clusters, ellipse = cluster_ellipse)

cluster_convex = [
    GO.convex_hull(Point.(tuple.(
        shapley_wide[shapley_wide.cluster .== cluster, :PC1], 
        shapley_wide[shapley_wide.cluster .== cluster, :PC2]
    )))
    for cluster in clusters
]

# fig_opts = (;
#     fontsize = fontsize,
#     size = (18.42centimetre, 18.42centimetre)
# )
# scale = scales(
#     X = (; label = "PC1"), 
#     Y = (; label = "PC2"),
#     Color = (; palette = [ESM_SSP_colors; vcat(guild_individual_colours...)])
# )
# legend_opts = (; position=:bottom, nbanks=3)

# ellipse = data(ESM_SSP_ellipse) * mapping(:ellipse, color=:ESM_SSP) * visual(Poly, alpha=0.3)
# lines = data(shapley_wide) * mapping(:PC1, :PC2, color=:output) * visual(Lines)

# fig = draw(ellipse + lines, scale; figure=fig_opts, legend=legend_opts)

fig = Figure(
    size = (28centimetre, 15centimetre),
    fontsize = fontsize
)
ax = Axis(fig[1,1], xlabel = "PC1", ylabel = "PC2", aspect = 1)

poly!(
    ax,
    ESM_SSP_convex,
    color=getindex.([Dict(ESM_SSP_colors)], ESM_SSP_ellipse.ESM_SSP),
    alpha=0.2
)
scatter!(
    ax, 
    shapley_wide.PC1, 
    shapley_wide.PC2, 
    color=getindex.([guild_individual_colours], shapley_wide.output),
    markersize=15,
    alpha=0.6
)

guild_elements = [MarkerElement(color=guild_individual_colours[g], marker=:circle) for g in guilds]
Legend(fig[2,1], guild_elements, getindex.([guild_clean_names], guilds), orientation=:horizontal, nbanks=5, title = "Guilds", tellwidth=false, colgap=8)

ESM_SSP_elements = [PolyElement(color=(Dict(ESM_SSP_colors)[esm_ssp], 0.6)) for esm_ssp in ESM_SSPs]
Legend(fig[3, 1], ESM_SSP_elements, getindex.([Dict(ESM_SSP_categories)], ESM_SSPs), title = "ESM - SSP", nbanks=2, patchsize = (10,5), orientation = :horizontal, tellwidth=false)
rowgap!(fig.layout, 2, Relative(0.005))


ax2 = Axis(fig[1,2], xlabel = "PC1", aspect = 1)
poly!(
    ax2,
    cluster_convex,
    color=getindex.([cluster_colours], cluster_ellipse.cluster),
    alpha=0.2
)
scatter!(
    ax2, 
    shapley_wide.PC1, 
    shapley_wide.PC2, 
    color=getindex.([guild_individual_colours], shapley_wide.output),
    markersize=15,
    alpha=0.6
)

cluster_elements = [PolyElement(color=(cluster_colours[cluster], 0.6)) for cluster in clusters]
Legend(fig[2, 2], cluster_elements, ["Boundary flow\ndriven", "Temperature and\nVertical mixing", "Vertical mixing\ndriven", "Temperature and\nboundary flows"], title = "Clusters", nbanks=2, patchsize = (10,5), orientation = :horizontal, tellwidth=false)
rowgap!(fig.layout, 2, Relative(0.005))


L = projection(M_pca)
ev = principalvars(M_pca) ./ sum(principalvars(M_pca))  # variance explained
corr_circle = L .* sqrt.(ev[1:2]')                      # scale loadings

scaling_factor = maximum(abs, y_pca)  # heuristic for visibility
L_scaled = corr_circle .* scaling_factor
corr_circle = L_scaled

ax3 = Axis(
    fig.layout[1,3],
    xticks = ax.xticks,
    xlabel = "PC1",
    # limits = (extrema(y_pca[1, :]), extrema(y_pca[2, :])),
    aspect=1
)
# Draw arrows for variables
# rename!(shapley_wide, "water_flows" => "boundary_flows")
for (i, variable) in enumerate(names(shapley_wide[:, Not(:output, :ESM_SSP, :cluster, :PC1, :PC2)]))
    arrows2d!(ax3, [mean(y_pca[1, :])], [mean(y_pca[2, :])], [corr_circle[i,1]], [corr_circle[i,2]], 
            shaftwidth=2, color=variable_colours[variable], alpha=0.6)
    # text!(ax, corr_circle[i,1], corr_circle[i,2], text=guilds[i], align=(:left, :bottom), alpha=0.5)
end

Label(fig.layout[1,1, TopLeft()], "A", font=:bold)
Label(fig.layout[1,2, TopLeft()], "B", font=:bold)
Label(fig.layout[1,3, TopLeft()], "C", font=:bold)

legend_entries = [PolyElement(color=variable_colours[variable]) for variable in variables]
Legend(
    fig[2,3],
    legend_entries,
    getindex.([variable_clean_names], variables),
    nbanks=2,
    patchsize = (7,10),
    label="Variable Groups",
    tellheight=false,
    tellwidth=false
)
rowgap!(fig.layout, 1, Relative(0.005))
linkyaxes!(ax, ax2, ax3)
linkxaxes!(ax, ax2, ax3)

save(joinpath(figs_path, "shapley_importance_clustered_pca.png"), fig, px_per_unit=dpi)