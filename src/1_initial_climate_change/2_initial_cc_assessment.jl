using CSV, DataFrames
using Parquet
using DataFrames
using Statistics
using MultivariateStats
using AlgebraOfGraphics
using GLMakie
using Colors

import GeometryBasics as GB

include("../analysis_common.jl")

CairoMakie.activate!()

master = CSV.read("../outputs/master_forcings_South_Africa_MA.csv", DataFrame)

result_files = readdir("../outputs/initial_runs/", join=true)
result_files = vcat(readdir.(result_files; join=true)...)
result_files = result_files[contains.(result_files, ["final_biomass"])]
result_files = result_files[contains.(result_files, "csv")]
variants = [first(match(r"(\d{4}-\d{4}-[[:upper:]]{4}-[[:lower:]]{3}\d{3})", fn).captures) for fn in result_files]

result_dfs = CSV.read.(result_files, [DataFrame])
for (r, result_df) in enumerate(result_dfs)
    result_df[!, :variant] .= variants[r]
end
result_df = vcat(result_dfs...)

# Plot timeseries of interesting biomass changes under climate change
result_df_lines = result_df[result_df.Description .∈ [guilds], :]
# result_df_lines = result_df_lines[result_df_lines.Description .∈ [result_df_lines.Description[contains.(result_df_lines.Description, ["fish"])]], :]
result_df_lines.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in result_df_lines.variant]
result_df_lines.ESM = [first(match(r"([[:upper:]]{4})", var).captures) for var in result_df_lines.variant]
result_df_lines.SSP = [first(match(r"([[:lower:]]{3}\d{3})", var).captures) for var in result_df_lines.variant]
result_df_lines.ESM_SSP = [first(match(r"([[:upper:]]{4}-[[:lower:]]{3}\d{3})", var).captures) for var in result_df_lines.variant]
result_df_lines.guild_clean = getindex.([guild_clean_names], result_df_lines.Description)

percent_change = combine(groupby(result_df_lines, [:Description, :ESM_SSP])) do sdf
    baseline = sdf[sdf.decade .== "2010-2019", :Model_annual_mean]
    (
        percent_change = (sdf.Model_annual_mean ./ baseline) .* 100,
        decade = sdf.decade
    )
end

fig_opts = (;
    fontsize = fontsize,
    size = (18.42centimetre, 18.42centimetre)
)
scale = scales(
    X = (; label = "Decade"), 
    Y = (; label = "Annual average biomass\n[mMN ⋅ gWW⁻¹]"),
    Color = (; label = "Earth System Model", categories = ["GFDL" => "GFDL-ESM4", "CNRM" => "CNRM-CM6-1-HR"]),
    LineStyle = (; label = "Socio-Economic Pathway", categories = ["ssp126" => "SSP1-2.6", "ssp370" => "SSP3-7.0"])
)
facet_opts = (; linkyaxes=:none)
legend_opts = (; position=:bottom)
axis_opts = (; xticklabelrotation = π/4)

biomass_timeseries = data(result_df_lines) * mapping(:decade, :Model_annual_mean, color=:ESM, linestyle=:SSP, layout=:guild_clean) * visual(Lines)
biomass_ts_fig = draw(biomass_timeseries, scale; facet=facet_opts, figure=fig_opts, legend=legend_opts, axis=axis_opts)

save("../figs/initial_cc_assessment/biomass_timeseries.png", biomass_ts_fig, px_per_unit=dpi)



# test = leftjoin(test, shapley_effects_wide[:, [:output, :cluster]], on = :Description => :output)
# test = sort(test, [:cluster, :Description])

# fig = Figure(size=(50centimetre, 20centimetre), fontsize=22pt)
# decades = unique(test.decade)
# guilds = unique(test.Description)
# clusters = unique(test[:, [:Description, :cluster]]).cluster
# guild_colours = Dict(zip(guilds, distinguishable_colors(length(guilds))))
# ax = Axis(fig[2,1], xlabel = "Decade", xticks = (1:length(decades), decades), ylabel = "Percentage of 2010-2019\n biomass [%]")

# test_data = test[test.ESM_SSP .== "CNRM-ssp370", :]
# guild_lines = [[(d, first(test_data[(test_data.Description .== guild) .& (test_data.decade .== decade), :].percent_change)) for (d, decade) in enumerate(decades)] for guild in guilds]
# map(x -> lines!(ax, guild_lines[x], color=cluster_colours[clusters[x]], linewidth=5), eachindex(guilds))
# Legend(fig[:,2], [PolyElement(color=cluster_colours[cluster]) for cluster in clusters], getindex.([guild_clean_names], guilds), "Guild Clusters")

# save("../figs/initial_cc_assessment/combined_cluster_importance_timeseries.png", fig, px_per_unit=dpi)

# PCA modelling of the biomass results
result_df_wide = unstack(result_df, :variant, :Description, :Model_annual_mean)
result_df_wide = result_df_wide[:, ["variant"; guilds]]
result_wide_standard = result_df_wide

rescaled_vars = guilds
for col in rescaled_vars
    if all(result_wide_standard[:, col] .== 0.0) continue end
    μ, σ = rescale!(result_wide_standard[!, col])
end

M_pca = fit(PCA, Matrix{Float64}(Matrix(result_df_wide[:, guilds])'); maxoutdim=2)
y_pca = predict(M_pca, Matrix{Float64}(Matrix(result_df_wide[:, guilds])'))
pca_val_df = DataFrame(variant = result_wide_standard.variant, PC1 = y_pca[1, :], PC2 = y_pca[2, :])

pca_val_df.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in pca_val_df.variant]
pca_val_df.ESM = [first(match(r"([[:upper:]]{4})", var).captures) for var in pca_val_df.variant]
pca_val_df.SSP = [first(match(r"([[:lower:]]{3}\d{3})", var).captures) for var in pca_val_df.variant]
pca_val_df.ESM_SSP = [first(match(r"([[:upper:]]{4}-[[:lower:]]{3}\d{3})", var).captures) for var in pca_val_df.variant]

fig_opts = (;
    fontsize = fontsize,
    size = (18.42centimetre, 14centimetre)
)
scale = scales(
    Color = (; label = "Earth System Model", categories = ["GFDL" => "GFDL-ESM4", "CNRM" => "CNRM-CM6-1-HR"]),
    LineStyle = (; label = "Socio-Economic Pathway", categories = ["ssp126" => "SSP1-2.6", "ssp370" => "SSP3-7.0"])
)
axis_opts = (; aspect=1)
legend_opts = (; position=:bottom, tellheight=false, tellwidth=false, nbanks=2)

# scat = data(pca_val_df) * mapping(:PC1, :PC2, marker=:decade, color=:ESM_SSP) * visual(Scatter, markersize=20)
line = data(pca_val_df) * mapping(:PC1, :PC2, color=:ESM, linestyle=:SSP) * visual(Lines)
pca_line_fig = draw(line, scale; figure=fig_opts, axis=axis_opts, legend=legend_opts)

fig = pca_line_fig.figure
ax = first(filter(x -> isa(x, Axis), fig.content))

arrowhead = GB.Polygon(Point2f[
    (0.15, 0),
    (-0.15,  0.12),
    (-0.15, -0.12)
])

for l in ax.scene.plots
    if l isa Makie.Lines
        points = l[1][]   # extract line data
        dx = points[end][1] - points[end-1][1]
        dy = points[end][2] - points[end-1][2]
        angle = atan(dy, dx)

        scatter!(
            ax,
            points[end],
            marker = :rtriangle,
            markersize = 15,
            rotation = angle,
            color = l.color[]
        )
        scatter!(
            ax,
            points[1],
            marker = :circle,
            markersize = 10,
            rotation = angle,
            color = l.color[]
        )
    end
end

L = projection(M_pca)
ev = principalvars(M_pca) ./ sum(principalvars(M_pca))  # variance explained
corr_circle = L .* sqrt.(ev[1:2]')                      # scale loadings

scaling_factor = maximum(abs, y_pca)  # heuristic for visibility
L_scaled = corr_circle .* scaling_factor
corr_circle = L_scaled

ax2 = Axis(
    fig.layout[1,2],
    xticks = ax.xticks,
    xlabel = "PC1",
    # limits = (extrema(y_pca[1, :]), extrema(y_pca[2, :])),
    aspect=1
)
# Draw arrows for variables
for (i, output) in enumerate(guilds)
    arrows2d!(ax2, [mean(y_pca[1, :])], [mean(y_pca[2, :])], [corr_circle[i,1]], [corr_circle[i,2]], 
            shaftwidth=2, color=guild_individual_colours[output], alpha=0.6)
    # text!(ax, corr_circle[i,1], corr_circle[i,2], text=guilds[i], align=(:left, :bottom), alpha=0.5)
end

Label(fig.layout[1,1, TopLeft()], "A", font=:bold)
Label(fig.layout[1,2, TopLeft()], "B", font=:bold)

legend_entries = [PolyElement(color=guild_individual_colours[guild]) for guild in guilds]
Legend(
    fig[2,2],
    legend_entries,
    guilds,
    nbanks=2,
    patchsize = (5,5),
    title="Guild",
    tellheight=false,
    tellwidth=false
)
rowsize!(fig.layout, 2, Relative(0.2))
rowgap!(fig.layout, 1, Relative(0.01))
linkyaxes!(ax, ax2)
linkxaxes!(ax, ax2)
fig

save("../figs/initial_cc_assessment/biomass_pca.png", fig, px_per_unit=dpi)

master.ESM_SSP = master.ESM .* "-" .* master.SSP
master_month_average = combine(groupby(master, [:variable, :decade, :SSP, :ESM]), :value => mean)

for vg in variables

    sub_drivers = master[master.variable .∈ [variable_groups[vg]], :]
    fig_opts = (;
        fontsize = fontsize,
        size = (18.42centimetre, 18.42centimetre)
    )
    # scale = scales(
    #     X = (; label = "Decade"), 
    #     Y = (; label = "Decadal Average"),
    #     Color = (; label = "Earth System Model", categories = ["GFDL" => "GFDL-ESM4", "CNRM" => "CNRM-CM6-1-HR"]),
    #     LineStyle = (; label = "Socio-Economic Pathway", categories = ["ssp126" => "SSP1-2.6", "ssp370" => "SSP3-7.0"])
    # )
    scale = scales(
        X = (; label = "Month"), 
        Y = (; label = "Monthly Value"),
        Color = (; label = "Decade"),
        LineStyle = (; label = "Earth System Model - Socio-Economic Pathway", categories = ["CNRM-ssp126" => "CNRM-CM6-1-HR\nSSP1-2.6", "CNRM-ssp370" => "CNRM-CM6-1-HR\nSSP3-7.0", "GFDL-ssp126" => "GFDL-ESM4\nSSP1-2.6", "GFDL-ssp370" => "GFDL-ESM4\nSSP3-7.0"])
    )
    facet_opts = (; linkyaxes=:none)
    legend_opts = (; position=:bottom, nbanks=3)
    axis_opts = (; xticklabelrotation = π/4)

    line = data(sub_drivers) * mapping(:month, :value, color=:decade, linestyle=:ESM_SSP, layout=:variable) * visual(Lines)
    fig = draw(line, scale; figure=fig_opts, axis=axis_opts, legend=legend_opts, facet=facet_opts)

    save("../figs/initial_cc_assessment/$(vg)_driver_timeseries.png", fig, px_per_unit=dpi)

end

# Plotting ecosystem mean trophic level over time
result_files = readdir("../outputs/initial_runs/", join=true)
result_files = vcat(readdir.(result_files; join=true)...)
result_files = result_files[contains.(result_files, ["ecosystem_indices"])]
result_files = result_files[contains.(result_files, "csv")]

indicators = ["wholeecosystemTL", "zooplanktonupTL", "fishupTL", "fishTL", "toppredTL"]
ecosystem_indices = [CSV.read(result_files[contains.(result_files, variant)], DataFrame) for variant in variants]
ecosystem_indices = [[parse(Float64, first(df[df.Description .== var, :].Value)) for var in indicators] for df in ecosystem_indices]

ecosystem_indices = DataFrame(hcat(variants, hcat(ecosystem_indices...)'), ["variant"; indicators])
ecosystem_indices = stack(ecosystem_indices, Not(:variant))

ecosystem_indices.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in ecosystem_indices.variant]
ecosystem_indices.ESM = [first(match(r"([[:upper:]]{4})", var).captures) for var in ecosystem_indices.variant]
ecosystem_indices.SSP = [first(match(r"([[:lower:]]{3}\d{3})", var).captures) for var in ecosystem_indices.variant]

fig_opts = (;
    fontsize = fontsize,
    size = (18.42centimetre, 18.42centimetre)
)
scale = scales(
    X = (; label = "Decade"), 
    Y = (; label = "Ecosystem Weighted Mean Trophic Level"),
    Color = (; label = "Earth System Model", categories = ["GFDL" => "GFDL-ESM4", "CNRM" => "CNRM-CM6-1-HR"]),
    LineStyle = (; label = "Socio-Economic Pathway", categories = ["ssp126" => "SSP1-2.6", "ssp370" => "SSP3-7.0"])
)
legend_opts = (; position=:bottom, orientation = :horizontal)
axis_opts = (; xticklabelrotation = π/4)
facet_opts = (; linkyaxes=:none)

line = data(ecosystem_indices) * mapping(:decade, :value, color=:ESM, linestyle=:SSP, layout=:variable) * visual(Lines)
fig = draw(line, scale; figure=fig_opts, axis=axis_opts, legend=legend_opts, facet=facet_opts)

save("../figs/initial_cc_assessment/ecosystem_mean_trophiclevel_timeseries.png", fig, px_per_unit=dpi)