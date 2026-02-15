"""
This script includes analyses and plotting of results from the base South_Africa_MA variants
(24 variants) capturing the impacts of climate change. This includes plotting basic timeseries
data, performing and plotting PCA, quantifying the degree of guild alignment with across ESM
and across SSP variation, plotting of various ecological indices, plotting diet compositions
and collating ecological information for demersal fish for discussions.
This script links with across_ESM_analysis and across_decade_analysis scripts as it defines
the esm_sep_guilds that are strongly associated with ESM variation and decade_sep_guilds
that are strongly associated with decade variation. These guilds are the targets for further
analyses.
"""

using CSV, DataFrames
using Parquet
using DataFrames
using Statistics
using MultivariateStats
using AlgebraOfGraphics
using GLMakie
using Colors
using LinearAlgebra

import GeometryBasics as GB

include("../analysis_common.jl")

CairoMakie.activate!()

monte_carlo_dir = "../outputs/initial_runs/monte_carlo_results/S_Benguela_MA/"

# 1. Plotting the initial assessment biomass timeseries
biomass_files = readdir("../outputs/initial_runs/", join=true)
biomass_files = biomass_files[isdir.(biomass_files)]
biomass_files = vcat(readdir.(biomass_files; join=true)...)
biomass_files = biomass_files[contains.(biomass_files, ["final_biomass"])]
biomass_files = biomass_files[contains.(biomass_files, "csv")]

variants = [first(match(r"(\d{4}-\d{4}-[[:upper:]]{4}-[[:lower:]]{3}\d{3})", fn).captures) for fn in biomass_files]

result_dfs = CSV.read.(biomass_files, [DataFrame])
for (r, result_df) in enumerate(result_dfs)
    result_df[!, :variant] .= variants[r]
end
result_df = vcat(result_dfs...)

#### Collate biomass data
result_df_lines = result_df[result_df.Description .∈ [guilds[guilds .!= "netprimprod"]], :]
result_df_lines.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in result_df_lines.variant]
result_df_lines.ESM = [first(match(r"([[:upper:]]{4})", var).captures) for var in result_df_lines.variant]
result_df_lines.SSP = [first(match(r"([[:lower:]]{3}\d{3})", var).captures) for var in result_df_lines.variant]
result_df_lines.ESM_SSP = [first(match(r"([[:upper:]]{4}-[[:lower:]]{3}\d{3})", var).captures) for var in result_df_lines.variant]
result_df_lines.guild_clean = getindex.([guild_clean_names], result_df_lines.Description)

result_df_lines[!, :min_quant] .= 0.0
result_df_lines[!, :max_quant] .= 0.0

result_df[!, :min_quant] .= 0.0
result_df[!, :max_quant] .= 0.0

for variant in variants
    mc_dir = joinpath(monte_carlo_dir, variant, "CredInt")
    aam_files = readdir(mc_dir; join=true)
    aam_file = first(aam_files[contains.(aam_files, ["AAMresults_whole-$(variant)-MC"])])

    aam_mc = CSV.read(aam_file, DataFrame)
    for guild in guilds
        if ismissing(se2e_mc_guilds[guild]) continue end
        
        guild_low_quant = first(aam_mc[aam_mc.Column1 .== "lowlimit", se2e_mc_guilds[guild]])
        result_df_lines[(result_df_lines.variant .== variant) .& (result_df_lines.Description .== guild), :min_quant] .= guild_low_quant
        result_df[(result_df.variant .== variant) .& (result_df.Description .== guild), :min_quant] .= guild_low_quant

        guild_upp_quant = first(aam_mc[aam_mc.Column1 .== "upplimit", se2e_mc_guilds[guild]])
        result_df_lines[(result_df_lines.variant .== variant) .& (result_df_lines.Description .== guild), :max_quant] .= guild_upp_quant
        result_df[(result_df.variant .== variant) .& (result_df.Description .== guild), :max_quant] .= guild_upp_quant
    end
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
biomass_bands_1 = data(result_df_lines[result_df_lines.SSP .== "ssp126", :]) * mapping(:decade, :min_quant, :max_quant, color=:ESM, layout=:guild_clean) * visual(Band; alpha = 0.2)
biomass_bands_2 = data(result_df_lines[result_df_lines.SSP .== "ssp370", :]) * mapping(:decade, :min_quant, :max_quant, color=:ESM, layout=:guild_clean) * visual(Band; alpha = 0.2)

biomass_ts_fig = draw(biomass_timeseries + biomass_bands_1 + biomass_bands_2, scale; facet=facet_opts, figure=fig_opts, legend=legend_opts, axis=axis_opts)
# biomass_ts_fig = draw(biomass_timeseries, scale; facet=facet_opts, figure=fig_opts, legend=legend_opts, axis=axis_opts)

save("../figs/initial_cc_assessment/biomass_timeseries_mc.png", biomass_ts_fig, px_per_unit=dpi)

# Additionally calculate percentage change data
percent_change = combine(groupby(result_df_lines, [:Description, :ESM_SSP, :ESM, :SSP, :guild_clean])) do sdf
    baseline_median = sdf[sdf.decade .== "2010-2019", :Model_annual_mean]
    baseline_min = sdf[sdf.decade .== "2010-2019", :min_quant]
    baseline_max = sdf[sdf.decade .== "2010-2019", :max_quant]
    (
        percent_change_median = ((sdf.Model_annual_mean .- baseline_median) ./ baseline_median ).* 100,
        percent_change_min = ((sdf.min_quant .- baseline_median) ./ baseline_median) .* 100,
        percent_change_max = ((sdf.max_quant .- baseline_median) ./ baseline_median) .* 100,
        decade = sdf.decade
    )
end

scale = scales(
    X = (; label = "Decade"), 
    Y = (; label = "Change in biomass from \n2010-2019 [%]"),
    Color = (; label = "Earth System Model", categories = ["GFDL" => "GFDL-ESM4", "CNRM" => "CNRM-CM6-1-HR"]),
    LineStyle = (; label = "Socio-Economic Pathway", categories = ["ssp126" => "SSP1-2.6", "ssp370" => "SSP3-7.0"])
)
percent_timeseries = data(percent_change) * mapping(:decade, :percent_change_median, color=:ESM, linestyle=:SSP, layout=:guild_clean) * visual(Lines)
percent_bands_1 = data(percent_change[percent_change.SSP .== "ssp126", :]) * mapping(:decade, :percent_change_min, :percent_change_max, color=:ESM, layout=:guild_clean) * visual(Band; alpha = 0.2)
percent_bands_2 = data(percent_change[percent_change.SSP .== "ssp370", :]) * mapping(:decade, :percent_change_min, :percent_change_max, color=:ESM, layout=:guild_clean) * visual(Band; alpha = 0.2)

percent_ts_fig = draw(percent_timeseries + percent_bands_1 + percent_bands_2, scale; figure=fig_opts, legend=legend_opts, axis=axis_opts)

save("../figs/initial_cc_assessment/percent_biomass_timeseries_mc.png", percent_ts_fig, px_per_unit=dpi)


# 2. plotting and analyses of biomass PCA
# Convert the data to a wide format for later ana
result_df_wide = unstack(result_df, :variant, :Description, :Model_annual_mean)
result_df_wide = result_df_wide[:, ["variant"; guilds[guilds .!= "netprimprod"]]]
result_wide_standard = result_df_wide
result_wide_standard.type .= "median"

result_df_min = unstack(result_df, :variant, :Description, :min_quant)
result_df_min = result_df_min[:, ["variant"; guilds[guilds .!= "netprimprod"]]]
result_df_min.type .= "min_quant"

result_df_max = unstack(result_df, :variant, :Description, :max_quant)
result_df_max = result_df_max[:, ["variant"; guilds[guilds .!= "netprimprod"]]]
result_df_max.type .= "max_quant"

result_wide_all = vcat(result_wide_standard, result_df_min, result_df_max)

# Rescale variables
rescaled_vars = guilds[guilds .!= "netprimprod"]
for col in rescaled_vars
    if all(result_wide_all[:, col] .== 0.0) continue end
    μ, σ = rescale!(result_wide_all[!, col])
end

# Perform PCA 
M_pca = fit(PCA, Matrix{Float64}(Matrix(result_wide_all[result_wide_all.type .== "median", rescaled_vars])'); maxoutdim=2)
y_pca = predict(M_pca, Matrix{Float64}(Matrix(result_wide_all[result_wide_all.type .== "median", rescaled_vars])'))

min_pca = predict(M_pca, Matrix{Float64}(Matrix(result_wide_all[result_wide_all.type .== "min_quant", rescaled_vars])'))
max_pca = predict(M_pca, Matrix{Float64}(Matrix(result_wide_all[result_wide_all.type .== "max_quant", rescaled_vars])'))

pca_val_df = DataFrame(
    variant = result_wide_standard.variant, 
    PC1 = y_pca[1, :], 
    PC2 = y_pca[2, :],
    min_PC1 = min_pca[1, :],
    min_PC2 = min_pca[2, :],
    max_PC1 = max_pca[1, :],
    max_PC2 = max_pca[2, :]
)

pca_val_df.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in pca_val_df.variant]
pca_val_df.ESM = [first(match(r"([[:upper:]]{4})", var).captures) for var in pca_val_df.variant]
pca_val_df.SSP = [first(match(r"([[:lower:]]{3}\d{3})", var).captures) for var in pca_val_df.variant]
pca_val_df.ESM_SSP = [first(match(r"([[:upper:]]{4}-[[:lower:]]{3}\d{3})", var).captures) for var in pca_val_df.variant]

# Quantify the separation between decades and between ESMs in PC space
esm_separation = (
    [mean(pca_val_df[pca_val_df.ESM .== "CNRM", :PC1]), mean(pca_val_df[pca_val_df.ESM .== "CNRM", :PC2])] .-
    [mean(pca_val_df[pca_val_df.ESM .== "GFDL", :PC1]), mean(pca_val_df[pca_val_df.ESM .== "GFDL", :PC2])]
)
decade_separation = (
    [mean(pca_val_df[pca_val_df.decade .== "2060-2069", :PC1]), mean(pca_val_df[pca_val_df.decade .== "2060-2069", :PC2])] .-
    [mean(pca_val_df[pca_val_df.decade .== "2050-2059", :PC1]), mean(pca_val_df[pca_val_df.decade .== "2050-2059", :PC2])] .-
    [mean(pca_val_df[pca_val_df.decade .== "2040-2049", :PC1]), mean(pca_val_df[pca_val_df.decade .== "2040-2049", :PC2])] .-
    [mean(pca_val_df[pca_val_df.decade .== "2030-2039", :PC1]), mean(pca_val_df[pca_val_df.decade .== "2030-2039", :PC2])] .-
    [mean(pca_val_df[pca_val_df.decade .== "2020-2029", :PC1]), mean(pca_val_df[pca_val_df.decade .== "2020-2029", :PC2])] .-
    [mean(pca_val_df[pca_val_df.decade .== "2010-2019", :PC1]), mean(pca_val_df[pca_val_df.decade .== "2010-2019", :PC2])]
)
pca_loadings = DataFrame(hcat(rescaled_vars, projection(M_pca)), ["guild", "PC1_loading", "PC2_loading"])

# Calculate the angle between PCA loadings and the separation vectors
function cos_similarity(a, b)
    return clamp(a⋅b/(norm(a)*norm(b)), -1, 1)
end

# Calculate angles to loadings
pca_loadings.esm_sep_similarity .= [abs(cos_similarity(esm_separation, [row.PC1_loading, row.PC2_loading])) for row in eachrow(pca_loadings)]
pca_loadings.decade_sep_similarity .= [abs(cos_similarity(decade_separation, [row.PC1_loading, row.PC2_loading])) for row in eachrow(pca_loadings)]
pca_loadings = sort(pca_loadings, :esm_sep_similarity, rev=true)
pca_loadings.guild_clean_names = getindex.([guild_clean_names], pca_loadings.guild)

# Define guilds that have an angle of less than 30 degrees
esm_sep_guilds = pca_loadings[pca_loadings.esm_sep_similarity .>= cosd(30), :guild]
decade_sep_guilds = pca_loadings[pca_loadings.decade_sep_similarity .>= cosd(30), :guild]

# Make plotting PCA confidence interval polygons
pca_polygons = combine(groupby(pca_val_df, [:ESM, :SSP])) do sdf
    min_points = GB.Point.(tuple.(sdf.min_PC1, sdf.min_PC2))
    max_points = GB.Point.(tuple.(sdf.max_PC1, sdf.max_PC2))
    
    all_points = vcat(min_points, max_points)
    hull = GO.convex_hull(all_points)
    (
        ;geometry = hull
    )
end

# Plot PCA data
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

line = data(pca_val_df) * mapping(:PC1, :PC2, color=:ESM, linestyle=:SSP) * visual(Lines)
pca_bands_1 = data(pca_polygons[pca_polygons.SSP .== "ssp126", :]) * mapping(:geometry, color=:ESM) * visual(Poly; alpha = 0.2)
pca_bands_2 = data(pca_polygons[pca_polygons.SSP .== "ssp370", :]) * mapping(:geometry, color=:ESM) * visual(Poly; alpha = 0.2)
pca_line_fig = draw(line + pca_bands_1 + pca_bands_2, scale; figure=fig_opts, axis=axis_opts, legend=legend_opts)

fig = pca_line_fig.figure
ax = first(filter(x -> isa(x, Axis), fig.content))

# Add arrowheads and points to the end and starts of the PCA lines to show time
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
        r_angle = atan(dy, dx)

        scatter!(
            ax,
            points[end],
            marker = :rtriangle,
            markersize = 15,
            rotation = r_angle,
            color = l.color[]
        )
        scatter!(
            ax,
            points[1],
            marker = :circle,
            markersize = 10,
            rotation = r_angle,
            color = l.color[]
        )
    end
end

# Plot PCA loading arrows
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
    aspect=1
)
# Draw arrows for variables
for (i, output) in enumerate(guilds[guilds .!= "netprimprod"])
    arrows2d!(ax2, [mean(y_pca[1, :])], [mean(y_pca[2, :])], [corr_circle[i,1]], [corr_circle[i,2]], 
            shaftwidth=2, color=guild_individual_colours[output], alpha=0.6)
end

lines!(ax2, [first(esm_separation), -first(esm_separation)], [last(esm_separation), -last(esm_separation)], linestyle=:dash, color=:red)
lines!(ax2, [first(decade_separation), -first(decade_separation)], [last(decade_separation), -last(decade_separation)], linestyle=:dash, color=:blue)

Label(fig.layout[1,1, TopLeft()], "A", font=:bold)
Label(fig.layout[1,2, TopLeft()], "B", font=:bold)

legend_entries = [PolyElement(color=guild_individual_colours[guild]) for guild in guilds[guilds .!= "netprimprod"]]
Legend(
    fig[2,2],
    legend_entries,
    guilds[guilds .!= "netprimprod"],
    nbanks=2,
    patchsize = (5,5),
    title="Guild",
    tellheight=false,
    tellwidth=false
)
Legend(
    fig[3,2], 
    [
        LineElement(color=:red, linestyle=:dash), 
        LineElement(color=:blue, linestyle=:dash)
    ],
    ["across ESM signal", "across decade signal"],
    tellheight=false,
    tellwidth=false,
    orientation=:horizontal
)
rowsize!(fig.layout, 2, Relative(0.2))
rowgap!(fig.layout, 1, Relative(0.01))
rowsize!(fig.layout, 3, Relative(0.05))
rowgap!(fig.layout, 2, Relative(0.05))
linkyaxes!(ax, ax2)
linkxaxes!(ax, ax2)

save("../figs/initial_cc_assessment/biomass_pca_mc.png", fig, px_per_unit=dpi)

# Plot the guild anlges to the separation vectors (small plot with thresholds indicated)
fig = Figure(
    fontsize = fontsize,
    size = (18.42centimetre, 10centimetre)
)
ax1 = Axis(fig[1,1], xlabel = "Contribution to ESM separation", ylabel = "Guilds", yticks=(1:nrow(pca_loadings), pca_loadings.guild_clean_names))
barplot!(ax1, 1:nrow(pca_loadings), pca_loadings.esm_sep_similarity; direction=:x)
vlines!(ax1, cosd(30), color=:red)

ax2 = Axis(fig[1,2], xlabel = "Contribution to decadal separation", yticksvisible=false, yticklabelsvisible=false)
barplot!(ax2, 1:nrow(pca_loadings), pca_loadings.decade_sep_similarity; direction=:x)
vlines!(ax2, cosd(30), color=:red)

save("../figs/initial_cc_assessment/decade_esm_separation_quantified.png", fig, px_per_unit=dpi)


# 3. Plotting ecosystem mean trophic level over time
indices_files = readdir("../outputs/initial_runs/", join=true)
indices_files = indices_files[isdir.(indices_files)]
indices_files = vcat(readdir.(indices_files; join=true)...)
indices_files = indices_files[contains.(indices_files, ["ecosystem_indices"])]
indices_files = indices_files[contains.(indices_files, "csv")]

indicators = ["wholeecosystemTL", "zooplanktonupTL", "fishupTL", "fishTL", "toppredTL"]
ecosystem_indices = [CSV.read(indices_files[contains.(indices_files, variant)], DataFrame) for variant in variants]
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

save("../figs/initial_cc_assessment/overall_trophiclevel_timeseries.png", fig, px_per_unit=dpi)


# 4. Plotting net primary production timeseries 
net_primprod = [CSV.read(indices_files[contains.(indices_files, variant)], DataFrame) for variant in variants]
net_primprod = [parse(Float64, first(df[df.Description .== "netprimprod", :].Value)) for df in net_primprod]

net_primprod = DataFrame(hcat(variants, hcat(net_primprod...)'), ["variant", "netprimprod"])

net_primprod.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in net_primprod.variant]
net_primprod.ESM = [first(match(r"([[:upper:]]{4})", var).captures) for var in net_primprod.variant]
net_primprod.SSP = [first(match(r"([[:lower:]]{3}\d{3})", var).captures) for var in net_primprod.variant]

fig_opts = (;
    fontsize = fontsize,
    size = (12centimetre, 10centimetre)
)
scale = scales(
    X = (; label = "Decade"), 
    Y = (; label = "Net Primary Production [mMN ⋅ m² ⋅ y]"),
    Color = (; label = "Earth System Model", categories = ["GFDL" => "GFDL-ESM4", "CNRM" => "CNRM-CM6-1-HR"]),
    LineStyle = (; label = "Socio-Economic Pathway", categories = ["ssp126" => "SSP1-2.6", "ssp370" => "SSP3-7.0"])
)
legend_opts = (; position=:bottom, orientation = :horizontal)

line = data(net_primprod) * mapping(:decade, :netprimprod, color=:ESM, linestyle=:SSP) * visual(Lines)
fig = draw(line, scale; figure=fig_opts, legend=legend_opts, facet=facet_opts)

save("../figs/initial_cc_assessment/net_primary_production_timeseries.png", fig, px_per_unit=dpi)


# 5. Plotting nitrogen flux data
flux_files = readdir("../outputs/initial_runs/", join=true)
flux_files = flux_files[.!contains.(flux_files, [".csv"])]
flux_files = vcat(readdir.(flux_files; join=true)...)
flux_files = flux_files[contains.(flux_files, ["flux_matrix_exclspawn"])]

flux_matrices = [CSV.read(flux_files[contains.(flux_files, variant)], DataFrame) for variant in variants]
for fmat in flux_matrices rename!(fmat, "Column1" => "source") end

# Extract total nitrogen flux before flux_matrices are later modified
total_flux = [DataFrame(guild = names(fmat[:, Not(:source)]),  total_uptake = vec(sum(Matrix(fmat[:, Not(:source)]); dims=1))) for fmat in flux_matrices]
for (v, df) in enumerate(total_flux) df.variant .= variants[v] end
total_flux = vcat(total_flux...)

# Proportion total nitrogen influx values by the guild annual mass
biomass = result_df_lines[result_df_lines.Description .∈ [decade_sep_guilds], [:Description, :variant, :Model_annual_mean]]
biomass.guild = getindex.([flux_guilds], biomass.Description)
total_flux = leftjoin(total_flux, biomass, on=[:guild, :variant])
total_flux.proportional_uptake = total_flux.total_uptake ./ total_flux.Model_annual_mean

# Take flux influx matrix values and convert them to proportions for each source guild (creating diet composition values)
function summarise_guild_fluxes(flux_matrix, sink_guild; guild_label=sink_guild)
    fluxes = flux_matrix[:, ["source", sink_guild]]
    total_influx = sum(fluxes[:, sink_guild])
    proportional_fluxes = fluxes[:, sink_guild] ./ total_influx

    result_df = DataFrame(source = fluxes.source, prop_influxes = proportional_fluxes)
    rename!(result_df, "prop_influxes" => "prop_influxes_$(guild_label)")

    return(result_df)
end

# summarise matrices across many sink guilds
function flux_matrix_extracted(flux_matrix, variants, f_index, guilds; flux_guilds=flux_guilds)
    guild_influxes = [summarise_guild_fluxes(flux_matrix, flux_guilds[g]; guild_label=g) for g in guilds]
    guild_influxes = reduce((df1, df2) -> leftjoin(df1, df2, on = :source), guild_influxes)
    guild_influxes = stack(guild_influxes; variable_name=:sink_guild, value_name=:influx)
    guild_influxes.variant .= variants[f_index]

    return guild_influxes
end

# Convert variant flux matrices to diet composition matrices
diet_matrices = [flux_matrix_extracted(fmat, variants, f, decade_sep_guilds) for (f, fmat) in enumerate(flux_matrices)]
diet_matrices = vcat(diet_matrices...)

## Plot the diet composition data as proportional bar plots
diet_comp = diet_matrices[diet_matrices.influx .!= 0.0, :]
diet_comp.sink_guild = last.(split.(diet_comp.sink_guild, ["prop_influxes_"]))
diet_comp.guild_clean = getindex.([guild_clean_names], diet_comp.sink_guild)
diet_comp.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in diet_comp.variant]
diet_comp.ESM = [first(match(r"([[:upper:]]{4})", var).captures) for var in diet_comp.variant]
diet_comp.SSP = [first(match(r"([[:lower:]]{3}\d{3})", var).captures) for var in diet_comp.variant]
diet_comp.ESM_SSP = diet_comp.ESM .* "-" .* diet_comp.SSP

source_guilds = unique(diet_comp.source)
source_guild_colours = [source_guilds[c] => col for (c, col) in enumerate(distinguishable_colors(length(source_guilds)))]

fig_opts = (; fontsize=fontsize, size=(18.42centimetre, 18.42centimetre))
scale = scales(
    X = (; label = "Sink guild"), 
    Y = (; label = "Proportion of nitrogen influx"),
    Color = (; label = "Source of nitrogen", palette = source_guild_colours),
    Col = (;categories = ESM_SSP_categories)
)
ax_opts = (; xticklabelrotation=π/2, ylabelpadding=10)

bars = data(diet_comp) * mapping(:guild_clean, :influx, row = :decade, col = :ESM_SSP, color=:source, stack=:source) * visual(BarPlot)
fig = draw(bars, scale; axis=ax_opts, figure=fig_opts)

save("../figs/initial_cc_assessment/decade_sep_guilds_nitrogen_sources.png", fig, px_per_unit=dpi)


# Plotting demersal fish and demersal fish larval influxes and outfluxes
# dfish_variant_biomass = result_df_lines[result_df_lines.Description .== "Demersal_fish_lar", [:variant, :Model_annual_mean]]
# rename!(dfish_variant_biomass, "Model_annual_mean" => "fy_biomass")

# dem_fish_influx = [df[:, ["source", "dfishlar"]] for df in flux_matrices]
# for (v, df) in enumerate(dem_fish_influx) df.variant .= variants[v] end
# dem_fish_influx = vcat(dem_fish_influx...)
# dem_fish_influx = leftjoin(dem_fish_influx, dfish_variant_biomass, on=:variant)
# dem_fish_influx.proportional_influx = dem_fish_influx.dfish ./ dem_fish_influx.fy_biomass

# dem_fish_outflux = [stack(df[df.source .== "dfish", :], Not(:source), variable_name=:sink, value_name=:outflux)[:, Not(:source)] for df in flux_matrices]
# for (v, df) in enumerate(dem_fish_outflux) df.variant .= variants[v] end
# dem_fish_outflux = vcat(dem_fish_outflux...)
# dem_fish_outflux = leftjoin(dem_fish_outflux, dfish_variant_biomass, on=:variant)
# dem_fish_outflux.proportional_outflux = dem_fish_outflux.outflux ./ dem_fish_outflux.fy_biomass


# plot_dfish_influx = dem_fish_influx[dem_fish_influx.dfishlar .> 0.0, :]
# plot_dfish_influx.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in plot_dfish_influx.variant]
# plot_dfish_influx.ESM = [first(match(r"([[:upper:]]{4})", var).captures) for var in plot_dfish_influx.variant]
# plot_dfish_influx.SSP = [first(match(r"([[:lower:]]{3}\d{3})", var).captures) for var in plot_dfish_influx.variant]
# fig_opts = (;
#     fontsize = fontsize,
#     size = (18.42centimetre, 14centimetre)
# )
# scale = scales(
#     Y = (; label = "Influx to Demersal fish from sources (proportional to dfish biomass)"),
#     Color = (; label = "Earth System Model", categories = ["GFDL" => "GFDL-ESM4", "CNRM" => "CNRM-CM6-1-HR"]),
#     LineStyle = (; label = "Socio-Economic Pathway", categories = ["ssp126" => "SSP1-2.6", "ssp370" => "SSP3-7.0"])
# )
# axis_opts = (; xticklabelrotation = π/4)
# legend_opts = (; position=:bottom, tellheight=false, tellwidth=false, nbanks=2)
# facet_opts = (; linkyaxes=false)

# # scat = data(pca_val_df) * mapping(:PC1, :PC2, marker=:decade, color=:ESM_SSP) * visual(Scatter, markersize=20)
# line = data(plot_dfish_influx) * mapping(:decade, :dfishlar, color=:ESM, linestyle=:SSP, layout=:source) * visual(Lines)
# fig = draw(line, scale; figure=fig_opts, axis=axis_opts, legend=legend_opts, facet=facet_opts)



# plot_dfish_outflux = dem_fish_outflux[dem_fish_outflux.outflux .> 0.0, :]
# plot_dfish_outflux.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in plot_dfish_outflux.variant]
# plot_dfish_outflux.ESM = [first(match(r"([[:upper:]]{4})", var).captures) for var in plot_dfish_outflux.variant]
# plot_dfish_outflux.SSP = [first(match(r"([[:lower:]]{3}\d{3})", var).captures) for var in plot_dfish_outflux.variant]
# fig_opts = (;
#     fontsize = fontsize,
#     size = (18.42centimetre, 14centimetre)
# )
# scale = scales(
#     Y = (; label = "Fluxes out of demersal fish to sinks (proportional to dfish biomass)"),
#     Color = (; label = "Earth System Model", categories = ["GFDL" => "GFDL-ESM4", "CNRM" => "CNRM-CM6-1-HR"]),
#     LineStyle = (; label = "Socio-Economic Pathway", categories = ["ssp126" => "SSP1-2.6", "ssp370" => "SSP3-7.0"])
# )
# axis_opts = (; xticklabelrotation = π/4)
# legend_opts = (; position=:bottom, tellheight=false, tellwidth=false, nbanks=2)
# facet_opts = (; linkyaxes=false)

# # scat = data(pca_val_df) * mapping(:PC1, :PC2, marker=:decade, color=:ESM_SSP) * visual(Scatter, markersize=20)
# line = data(plot_dfish_outflux) * mapping(:decade, :proportional_outflux, color=:ESM, linestyle=:SSP, layout=:sink) * visual(Lines)
# fig = draw(line, scale; figure=fig_opts, axis=axis_opts, legend=legend_opts, facet=facet_opts)



# fig_opts = (;
#     fontsize = fontsize,
#     size = (18.42centimetre, 14centimetre)
# )
# scale = scales(
#     Y = (; label = "Total nitrogen uptake relative to guild biomass"),
#     Color = (; label = "Earth System Model", categories = ["GFDL" => "GFDL-ESM4", "CNRM" => "CNRM-CM6-1-HR"]),
#     LineStyle = (; label = "Socio-Economic Pathway", categories = ["ssp126" => "SSP1-2.6", "ssp370" => "SSP3-7.0"])
# )
# axis_opts = (; xticklabelrotation = π/4)
# legend_opts = (; position=:bottom, tellheight=false, tellwidth=false, nbanks=2)
# facet_opts = (; linkyaxes=false)

# # scat = data(pca_val_df) * mapping(:PC1, :PC2, marker=:decade, color=:ESM_SSP) * visual(Scatter, markersize=20)
# line = data(total_flux) * mapping(:decade, :proportional_uptake, color=:ESM, linestyle=:SSP, layout=:guild) * visual(Lines)
# fig = draw(line, scale; figure=fig_opts, axis=axis_opts, legend=legend_opts, facet=facet_opts)


# 6. Collate the trophic production data (biological rates) for the demersal fish and demersal fish larvae guilds
# This information is saved out to csvs for creating tables for communication of high turnover rates of dfish

prod_files = readdir("../outputs/initial_runs/", join=true)
prod_files = prod_files[isdir.(prod_files)]
prod_files = vcat(readdir.(prod_files; join=true)...)
prod_files = prod_files[contains.(prod_files, ["trophic_prod"])]

prod_files = [CSV.read(prod_files[contains.(prod_files, variant)], DataFrame) for variant in variants]
prod_files = [df[:, [:Column1, :dfish, :dfishlar]] for df in prod_files]
for (v, df) in enumerate(prod_files) rename!(df, "dfish" => "Demersal fish $(variants[v])", "dfishlar" => "Demersal fish larvae $(variants[v])") end

wide_trophic_prod = reduce((df1, df2) -> leftjoin(df1, df2, on = :Column1), prod_files)
rename!(wide_trophic_prod, "Column1" => "Rate")
wide_trophic_prod = wide_trophic_prod[wide_trophic_prod.Rate .∈ [[
        "T_uptake", 
        "T_defec_excr", 
        "Recruitment", 
        "Spawning", 
        "Net_Prod_feeding", 
        "Net_Prod"
    ]],
:]

CSV.write("../outputs/initial_runs/demersal_fish_rates_variants.csv", wide_trophic_prod)

desired_col_labels = ["Demersal fish 2010-2019-CNRM-ssp370",
    "Demersal fish larvae 2010-2019-CNRM-ssp370",
    "Demersal fish 2010-2019-GFDL-ssp370",
    "Demersal fish larvae 2010-2019-GFDL-ssp370"
]
sub_wide_dfish = wide_trophic_prod[:, ["Rate"; desired_col_labels]]
for col in desired_col_labels
    sub_wide_dfish[!, col] = round.(sub_wide_dfish[:, col], sigdigits=4)
end
CSV.write("../outputs/initial_runs/demersal_fish_rates_variants_filtered.csv", sub_wide_dfish)


# Demersal fish parameter sensitivity
target_guilds = ["Planktivorous_fish", "Demersal_fish", "Demersal_fish_larvae", "Birds", "Cetaceans"]

param_permutations = CSV.read("../outputs/initial_runs/demfish_sens_params.csv", DataFrame)[:, 2:5]
param_permutations = stack(param_permutations, Not([:param_name, :param_id]), variable_name = :ESM, value_name = :param_value)
sens_variants = ["2010-2019-CNRM-ssp370", "2010-2019-GFDL-ssp370"]

base_params = CSV.read("../outputs/initial_runs/demfish_sens_base_params.csv", DataFrame)
base_params = stack(base_params, Not(:param_name), variable_name = :ESM, value_name = :param_value)

for tg in target_guilds
    param_permutations[!, tg] .= 0.0
    base_params[!, tg] .= 0.0
end

for variant in sens_variants
    esm = first(match(r"([[:upper:]]{4})", variant).captures)
    perm_files = readdir("../outputs/initial_runs/$(variant)/demfish_sens/"; join = true)
    perm_files = perm_files[contains.(perm_files, ["final_year"])]

    for pn in unique(param_permutations.param_id)
        pfn = first(perm_files[contains.(perm_files, ["p$(pn).csv"])])
        final_year = CSV.read(pfn, DataFrame)
        
        for tg in target_guilds
            tg_mass = first(final_year[final_year.Description .== tg, :].Model_annual_mean)
            param_permutations[(param_permutations.ESM .== esm) .& (param_permutations.param_id .== pn), tg] .= tg_mass
        end
    end

    for tg in target_guilds
        tg_mean_mass = first(result_df_lines[(result_df_lines.Description .== tg) .& (result_df_lines.variant .== variant), :].Model_annual_mean)
        base_params[base_params.ESM .== esm, tg] .= tg_mean_mass
    end
end

param_permutations = stack(param_permutations, Not([:param_name, :ESM, :param_value, :param_id]), variable_name = :guild, value_name = :biomass)
param_permutations.guild_clean_name = getindex.([guild_clean_names], param_permutations.guild)
base_params = stack(base_params, Not([:param_name, :ESM, :param_value]), variable_name = :guild, value_name = :biomass)
base_params.guild_clean_name = getindex.([guild_clean_names], base_params.guild)

fig_opts = (;
    fontsize = fontsize,
    size = (18.42centimetre, 18.42centimetre)
)
facet_opts = (; linkxaxes=:none)
legend_opts = (; position=:bottom)
axis_opts = (; xticklabelrotation = π/4)

# Plots for larval parameter sensitivity analysis
lar_params_scale = scales(
    X = (; label = "Parameter value"), 
    Y = (; label = "Annual average biomass\n[mMN ⋅ gWW⁻¹]"),
    Color = (; label = "Earth System Model", categories = ESM_categories),
    Row = (; categories = param_categories[contains.(first.(param_categories), ["lar"])])
)

demfishlar_sens_scatter = data(param_permutations[contains.(param_permutations.param_name, ["lar"]), :]) * mapping(:param_value, :biomass, color=:ESM, row=:param_name, col=:guild_clean_name) * visual(Scatter; alpha=0.4)
demfishlar_sens_means = data(base_params[contains.(base_params.param_name, ["lar"]), :]) * mapping(:param_value, :biomass, row=:param_name, col = :guild_clean_name) * visual(Scatter, alpha = 0.4, marker=:utriangle, color=:red, markersize=15)

demfishlar_sens_fig = draw(demfishlar_sens_scatter + demfishlar_sens_means, lar_params_scale; facet=facet_opts, figure=fig_opts, legend=legend_opts, axis=axis_opts)
save("../figs/initial_cc_assessment/demfishlar_param_sens.png", demfishlar_sens_fig, px_per_unit=dpi)

# Plots for adult parameter sensitivity analysis
dfish_params_scale = scales(
    X = (; label = "Parameter value"), 
    Y = (; label = "Annual average biomass\n[mMN ⋅ gWW⁻¹]"),
    Color = (; label = "Earth System Model", categories = ESM_categories),
    Row = (; categories = param_categories[.!contains.(first.(param_categories), ["lar"])])
)

dfish_sens_scatter = data(param_permutations[.!contains.(param_permutations.param_name, ["lar"]), :]) * mapping(:param_value, :biomass, color=:ESM, row=:param_name, col=:guild_clean_name) * visual(Scatter; alpha=0.4)
dfish_sens_means = data(base_params[.!contains.(base_params.param_name, ["lar"]), :]) * mapping(:param_value, :biomass, row=:param_name, col = :guild_clean_name) * visual(Scatter, alpha = 0.4, marker=:utriangle, color=:red, markersize=15)

dfish_sens_fig = draw(dfish_sens_scatter + dfish_sens_means, dfish_params_scale; facet=facet_opts, figure=fig_opts, legend=legend_opts, axis=axis_opts)
save("../figs/initial_cc_assessment/demfish_param_sens.png", dfish_sens_fig, px_per_unit=dpi)