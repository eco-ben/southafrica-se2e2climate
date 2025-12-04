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


# 1. Plotting Environmental driving data timeseries
master = CSV.read("../outputs/master_forcings_South_Africa_MA.csv", DataFrame)

master.ESM_SSP = master.ESM .* "-" .* master.SSP
master_month_average = combine(groupby(master, [:variable, :decade, :SSP, :ESM]), :value => mean)

for vg in variables
    sub_drivers = master[master.variable .∈ [variable_groups[vg]], :]
    fig_opts = (;
        fontsize = fontsize,
        size = (18.42centimetre, 18.42centimetre)
    )
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

    save("../figs/env_drivers/$(vg)_driver_timeseries.png", fig, px_per_unit=dpi)
end

# 2. Plotting the average variance of each driving data variable across decades and ESMs analysed.
function grouped_coefficient_of_var(values, groupings)
    Dict(k => Float64[] for k in unique(groupings))
end

annual_variance = combine(groupby(master, [:variable, :decade, :SSP, :ESM]), :value => var)
annual_variance = annual_variance[annual_variance.value_var .!= 0.0, :]

across_decade = combine(groupby(annual_variance, [:variable, :decade]), :value_var => mean)
across_decade_cv = combine(groupby(across_decade, :variable), :value_var_mean => variation)
across_decade_cv[!, "analysis_level"] .= "across decade"

across_esm = combine(groupby(annual_variance, [:variable, :ESM]), :value_var => mean)
across_esm_cv = combine(groupby(across_esm, :variable), :value_var_mean => variation)
across_esm_cv[!, "analysis_level"] .= "across esm"

coef_variation_data = vcat(across_decade_cv, across_esm_cv)
coef_variation_data = coef_variation_data[coef_variation_data.variable .∈ [vcat(collect(values(variable_groups))...)], :]
coef_variation_data.variable_group = [first([k for (k, v) in variable_groups if in(var, v)]) for var in coef_variation_data.variable]
coef_variation_data.variable_clean_group = getindex.([variable_clean_names], coef_variation_data.variable_group)

fig_opts = (;
    fontsize = fontsize,
    size = (18centimetre, 14centimetre)
)
scale = scales(
    X = (; label = "Individual variables"), 
    Y = (; label = "Coefficient of Variation across levels"),
    Color = (; label = "Variable groups"),
    Layout = (; categories = ["across decade" => "Across decades", "across esm" => "Across ESMs"])
)
axis_opts = (; xticklabelrotation = π/2)
legend_opts = (; position=:bottom, orientation = :horizontal, nbanks=2)

bar = data(sort(coef_variation_data, :variable_group)) * mapping(:variable => presorted, :value_var_mean_variation , color=:variable_clean_group, layout=:analysis_level) * visual(BarPlot)
fig = draw(bar, scale; figure = fig_opts, legend = legend_opts, axis = axis_opts)

save("../figs/env_drivers/driving_variable_variation.png", fig, px_per_unit=dpi)
