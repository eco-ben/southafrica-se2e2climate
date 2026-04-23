using CSV, DataFrames
using Statistics
using MultivariateStats
using AlgebraOfGraphics
using CairoMakie
using Clustering
using CategoricalArrays
using Distances

include("../analysis_common.jl")

output_path = "../outputs/across_decade_permutations"
figs_path = "../figs/across_decade_permutations"

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

function collate_variable_interactions(permutation_plan, perm_plan_outputs, base_variant, comp_variant, variable, variables, guilds; permutation_id_col=:perm_id)
    others = setdiff(variables, [variable])
    interactions = vcat.(variable, others)
    
    variable_comparison_levels = DataFrame(variable_changed = vcat("baseline", [[variable]], interactions))
    variable_comparison_levels.perm_id = map(
        x -> find_permutation(permutation_plan, base_variant, comp_variant, x, variables),
        variable_comparison_levels.variable_changed
    )

    for guild in guilds
        variable_comparison_levels[!, guild] = perm_plan_outputs[variable_comparison_levels.perm_id, guild]
        # variable_comparison_levels[!, "$(guild)_percent_change"] = 100 .- (variable_comparison_levels[variable_comparison_levels.variable_changed .== "baseline", guild] ./ 
        #     variable_comparison_levels[:, guild] .* 100)
        variable_comparison_levels[!, "$(guild)_percent_change"] = (variable_comparison_levels[:, guild] .- variable_comparison_levels[variable_comparison_levels.variable_changed .== "baseline", guild]) ./ 
            variable_comparison_levels[variable_comparison_levels.variable_changed .== "baseline", guild] .* 100
    end

    variable_comparison_levels.labels = ifelse.(
        variable_comparison_levels.variable_changed .== "baseline",
        variable_comparison_levels.variable_changed,
        join.(variable_comparison_levels.variable_changed, ":")
    )
    variable_comparison_levels.analysis_variable .= variable

    return variable_comparison_levels
end

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

# 1. Plot all Shapley Effects for all guilds and all ESM-SSP combinations
shapley_effects_all = [CSV.read(joinpath(output_path, esm_ssp, "shapley_effects.csv"), DataFrame) for esm_ssp in ESM_SSPs]
for (e, df) in enumerate(shapley_effects_all) df[!, "ESM_SSP"] .= ESM_SSPs[e] end
shapley_effects_all = vcat(shapley_effects_all...)


# 2. Plotting Shapley Effect contributions for guilds that contribute most to across decade signals
# This includes acrosss decade Shapley Effect values for all ESM-SSP combinations

all_shap_effect = shapley_effects_all
rename!(all_shap_effect, ["output" => "guild", "variable_group" => "variable"])
all_shap_effect.variable = ifelse.(all_shap_effect.variable .== "water_flows", "boundary_flows", all_shap_effect.variable)
all_shap_effect = all_shap_effect[all_shap_effect.variable .!= "constant", :]

decade_sep_shap = all_shap_effect[all_shap_effect.guild .∈ [["netprimprod"; decade_sep_guilds]], :]
decade_sep_shap = sort(decade_sep_shap, :variable, rev=true)

guild_order = Dict(
    "netprimprod" => 1,
    "Deep_layer_phytoplankton" => 2,
    "Benthos_carn/scav_feeders" => 3,
    "Demersal_fish" => 4,
    "Pinnipeds" => 5,
    "Birds" => 6,
    "Cetaceans" => 7
)
decade_sep_shap = sort(decade_sep_shap, :guild, by = x -> guild_order[x])

decade_sep_shap.variable_clean_name = getindex.([variable_clean_names], decade_sep_shap.variable)
decade_sep_shap.guild_clean_name = getindex.([guild_clean_names], decade_sep_shap.guild)
decade_sep_shap.ESM = first.(split.(decade_sep_shap.ESM_SSP, ["-"]))
decade_sep_shap.SSP = last.(split.(decade_sep_shap.ESM_SSP, ["-"]))

variable_colours = [
    "vertical mixing" => Makie.wong_colors()[4],
    "temperature" => Makie.wong_colors()[6],
    "river outputs" => :grey,
    "nutrient concentrations" => :grey,
    "light" => :grey,
    "boundary flows" => Makie.wong_colors()[1],
    "atmospheric nutrient flux" => :grey
]

fig_opts = (; fontsize=fontsize, size=(18.42centimetre, 18.42centimetre))
scale = scales(
    X = (; label = " "), 
    Y = (; label = "Shapley Effect"),
    Color = (; label = "Variable group", palette = variable_colours),
    Row = (; categories = ESM_categories),
    Col = (; categories = SSP_categories)
)
ax_opts = (; xticklabelrotation=π/4, ylabelpadding=10)

bars = data(decade_sep_shap) * mapping(:guild_clean_name => sorter(decade_sep_shap.guild_clean_name), :shapley_effect, row=:ESM, col=:SSP, color=:variable_clean_name, stack=:variable) * visual(BarPlot)
fig = draw(bars, scale; axis=ax_opts, figure=fig_opts)

save("../figs/across_decade_permutations/important_guild_shapley.png", fig, px_per_unit=dpi)


# 3. Plotting variable main and interaction effects using across decade permutations
"""
    get_esmssp_interactions(esm_ssp, output_path; baseline_decade="2010-2019", comparison_decade="2060-2069")

For a single ESM-SSP combination, using the permutation lookup table this function extracts
the guild outputs for each variable and calculates the effect size for main variable effects
and variable interaction effects.
Effect sizes are calculated as the percentage change in biomass compared to a baseline value.
The baseline value is the permutation where all variables are at the baseline decade.
Main effects are calculated where a single variable is changed from the baseline to the comparison
decade.
Variable interaction effects are calculated where a pair of variables are changed from the
baseline to the comparison decade.

# Arguments
- `esm_ssp` : ESM-SSP combination for analysis. The function can only run on a single ESM-SSP combination.
- `output_path` : base output path for across decade permutations.
- `baseline_decade` : baseline decade for comparisons, default value is 2010-2019.
- `comparison_decade` : default value is 2060-2069. If it is the same as the baseline decade all effects will be 0.

# Returns
Tuple containing a dataframe with the main effects and a dataframe with the variable interaction effects.
"""
function get_esmssp_interactions(esm_ssp, output_path; baseline_decade="2010-2019", comparison_decade="2060-2069")
    sub_output_path = joinpath(output_path, esm_ssp)

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
        if any(perm_outputs[1].Description .== guild)
            output = [out_df[out_df.Description .== [guild], :Model_annual_mean] for out_df in perm_outputs]
            output = vcat(output...)
            perm_plan_outputs[!, guild] = output
        end
    end

    files = readdir(sub_output_path; join=true)
    files = files[contains.(files, ["ecosystem_indices_"])]
    perm_outputs = [CSV.read(files[contains.(files, ["perm_$(pid).csv"])], DataFrame) for pid in permutation_plan.perm_id]
    
    net_prim_prod = vcat([parse.(Float64, out_df[out_df.Description .== ["netprimprod"], :Value]) for out_df in perm_outputs]...)
    perm_plan_outputs[!, "netprimprod"] = net_prim_prod
    
    base_variant = baseline_decade * "-" * esm_ssp
    comp_variant = comparison_decade * "-" * esm_ssp

    var_permutation_outputs = vcat([
        collate_variable_interactions(permutation_plan, perm_plan_outputs, base_variant, comp_variant, var, variables, guilds)
        for var in variables
    ]...)
    percentage_columns = guilds .* "_percent_change"
    
    var_permutation_outputs = stack(var_permutation_outputs, Not(:variable_changed, :perm_id, :labels, :analysis_variable), variable_name=:guild, value_name=:value)
    var_permutation_outputs = var_permutation_outputs[var_permutation_outputs.variable_changed .!= "baseline", :]
    
    base_variable_changes = var_permutation_outputs[.!occursin.(":", var_permutation_outputs.labels), :]
    interaction_rows = var_permutation_outputs[
        occursin.(":", var_permutation_outputs.labels) .& occursin.("percent_change", var_permutation_outputs.guild), 
    :]
    
    variable_interaction_labels = combine(groupby(interaction_rows, [:labels, :analysis_variable, :guild])) do sdf
        dfrow = sdf[1, :]

        base_variable = first(split(dfrow.labels, ":"))
        interacting_variable = last(split(dfrow.labels, ":"))
        guild = dfrow.guild

        Δ_base_var = first(base_variable_changes[
            (base_variable_changes.labels .== base_variable) .& 
            (base_variable_changes.guild .== guild), 
        :value])
        Δ_second_var = first(base_variable_changes[
            (base_variable_changes.labels .== interacting_variable) .& 
            (base_variable_changes.guild .== guild), 
        :value])

        direction_type = sign(Δ_base_var) == sign( Δ_second_var) ? "same_direction" : "opposing_direction"
        
        expected_additive_change = Δ_base_var + Δ_second_var
        tolerance = max(0.005, 0.05 * abs(expected_additive_change))
        Δ_interaction = dfrow.value
        residual = Δ_interaction - expected_additive_change

        interaction_type = abs(residual) <= tolerance ? "linear" : "non-linear"

        (;
            direction_class = direction_type,
            interaction_class = interaction_type,
            interaction_value = Δ_interaction,
            pure_interaction = residual,
            interaction_ratio = abs(Δ_interaction) != 0 ? abs(residual) / abs(expected_additive_change) : 0
        )
    end
    
    base_variable_changes.ESM_SSP .= esm_ssp
    variable_interaction_labels.ESM_SSP .= esm_ssp

    return (main_effects = base_variable_changes, interaction_effects = variable_interaction_labels)
end

var_main_effects = vcat(first.(get_esmssp_interactions.(ESM_SSPs, output_path))...) 
var_interaction_data = vcat(last.(get_esmssp_interactions.(ESM_SSPs, output_path))...)

normalize_label(s) = join(sort(split(s, ":")), ":")
filtered = combine(groupby(var_interaction_data, [:ESM_SSP, :guild])) do sdf
    df = DataFrames.transform(sub, :labels => ByRow(normalize_label) => :label_normalized)
    df = unique(df[:, [:guild, :ESM_SSP, :label_normalized, :pure_interaction]])

    (label = df.label_normalized, interaction = df.pure_interaction)
end

# Plotting the distrubution of variable interaction effects for all guilds
fig = Figure(size = (10centimetre, 10centimetre), fontsize=fontsize)
ax = Axis(fig[1,1], xlabel = "Environmental variable group interaction effect \n[% change in biomass from 2010-2019 baseline]", ylabel="number of samples")
hist!(ax, filtered.interaction, bins=50, color = :red, strokewidth = 1, strokecolor = :black)
vlines!(ax, quantile.([filtered.interaction], [0.1, 0.9]), color=:blue)
save("../figs/across_decade_permutations/variable_interaction_histogram.png", fig, px_per_unit=dpi)