using CSV, DataFrames
using Parquet
using DataFrames
using Statistics
using MultivariateStats
using AlgebraOfigraphics
using GLMakie
using Colors
using StatsBase

import GeometryBasics as GB

include("analysis_common.jl")

CairoMakie.activate!()


# 1. Plotting Environmental driving data timeseries
master = CSV.read("../outputs/master_forcings_South_Africa_MA.csv", DataFrame)

master.ESM_SSP = master.ESM .* "-" .* master.SSP
master_month_average = combine(groupby(master, [:variable, :decade, :SSP, :ESM]), :value => mean)

for vg in variables
    sub_drivers = master[master.variable .∈ [variable_groups[vg]], :]
    sub_drivers = sub_drivers[sub_drivers.decade .!= "2010-2015", :]
    sub_drivers.variable_label = getindex.([variable_units], sub_drivers.variable)
    sub_drivers.variable_label = sub_drivers.variable .* " [" .* sub_drivers.variable_label .* "]"

    fig_opts = (;
        fontsize = fontsize,
        size = (18.42centimetre, 18.42centimetre)
    )
    scale = scales(
        X = (; label = "Month"), 
        Y = (; label = "Monthly Value"),
        LineStyle = (; label = "SSP", categories = SSP_categories, palette=SSP_linestyles),
        Color = (; label = "NEMO-ERSEM forcing model", categories = ESM_categories)
    )
    facet_opts = (; linkyaxes=:minimal)
    legend_opts = (; position=:bottom, nbanks=3)

    axis_opts = (; xticks=([1, 6, 12], ["Jan", "Jun", "Dec"]), yticks=LinearTicks(3))

    line = data(sub_drivers) * mapping(:month, :value, color=:ESM, linestyle=:SSP, col=:decade, row=:variable_label) * visual(Lines)
    fig = draw(line, scale; figure=fig_opts, axis=axis_opts, legend=legend_opts, facet=facet_opts)

    if vg ∈ ["nutrient_concentrations", "atmospheric_nutrient_flux"]
        right_col = size(fig.figure.layout)[2]

        for gc in fig.figure.layout.content
            # Check it's a Label and that it lives in the rightmost column (row labels)
            if (gc.content isa Makie.Label) && (gc.span.cols == right_col:right_col) && (gc.content.rotation[] ≈ -π/2)
                gc.content.rotation[] = 0   # 45 degrees counter-clockwise
            end
        end
    end
    

    save("../figs/env_drivers/$(vg)_driver_timeseries.png", fig, px_per_unit=dpi)
end

# 2. Plotting the average variance of each driving data variable across decades and ESMs analysed.

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
    Layout = (; categories = ["across decade" => "Across decades", "across esm" => "Across NEMO-ERSEM forcing models"])
)
axis_opts = (; xticklabelrotation = π/2)
legend_opts = (; position=:bottom, orientation = :horizontal, nbanks=2)

bar = data(sort(coef_variation_data, :variable_group)) * mapping(:variable => presorted, :value_var_mean_variation , color=:variable_clean_group, layout=:analysis_level) * visual(BarPlot)
fig = draw(bar, scale; figure = fig_opts, legend = legend_opts, axis = axis_opts)

save("../figs/env_drivers/driving_variable_variation.png", fig, px_per_unit=dpi)
