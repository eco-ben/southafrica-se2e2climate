using GlobalSensitivity
using CSV, DataFrames
using Parquet
using DataFrames
using Combinatorics
using Statistics
using MultivariateStats
using MLDataUtils

master = CSV.read("../outputs/master_forcings_South_Africa_MA.csv", DataFrame)

function find_encoded_val(full_code_value, encoding_scheme, disregarded_pattern)
    split_code_value = replace(full_code_value, disregarded_pattern => "")

    return encoding_scheme[split_code_value]
end
function extract_annual_results(filename, id_column, id, value_column)
    df = DataFrame(read_parquet(filename))

    return first(df[df[:, id_column] .== id, value_column])
end

within_decade_perms = CSV.read("../outputs/within_decade_ssp_esm_permutations.csv", DataFrame)
input_variables = names(within_decade_perms)[(names(within_decade_perms) .!= "perm_id") .& (names(within_decade_perms) .!= "Column1")]
ssp_esm_encoding = Dict("CNRM-ssp126" => 1, "GFDL-ssp126" => 2, "CNRM-ssp370" => 3, "GFDL-ssp370" => 4)

# For 2010-2019:
decade_perms_2010_2019 = within_decade_perms[contains.(within_decade_perms.light, "2010-2019"), :]
for var in input_variables
    decade_perms_2010_2019[!, var] = [find_encoded_val(x, ssp_esm_encoding, r"(\d{4}-\d{4}-)") for x in decade_perms_2010_2019[:, var]]
end

result_files = 1:200
result_files = ["../outputs/within_decade_permutations/model_outputs_perm_$(x).parq" for x in result_files]

df = decade_perms_2010_2019[1:200, :]
df.annual_surface_phyt = [extract_annual_results(filename, "Description", "Surface_layer_phytoplankton", "Model_annual_mean") for filename in result_files]
df.annual_dem_fish = [extract_annual_results(filename, "Description", "Demersal_fish", "Model_annual_mean") for filename in result_files]
df.annual_plank_fish = [extract_annual_results(filename, "Description", "Planktivorous_fish", "Model_annual_mean") for filename in result_files]
df.annual_nit_mass = [extract_annual_results(filename, "Description", "Total_nitrogen_mass", "Model_annual_mean") for filename in result_files]
df.annual_birds = [extract_annual_results(filename, "Description", "Birds", "Model_annual_mean") for filename in result_files]
outputs = [
    :annual_surface_phyt, 
    :annual_dem_fish, 
    :annual_plank_fish, 
    :annual_nit_mass,
    :annual_birds
]

function variable_contribution_weighted(df, variable, output_col)
    groups = groupby(df, variable)
    group_sizes = [nrow(g) for g in groups]
    group_means = [mean(g[:, output_col]) for g in groups]
    
    total_mean = mean(df[:, output_col])
    
    # Weighted variance
    variable_variance = sum(w * (m - total_mean)^2 for (w,m) in zip(group_sizes, group_means)) / sum(group_sizes)
    
    # Normalize
    scaled_variable_variance = variable_variance / Statistics.var(df[:, output_col]; corrected=false)
    return scaled_variable_variance
end

"""
    v_of_S(df::DataFrame, S::Union{Vector{Symbol}, Vector{String}}, y::Symbol)

Calculate the variance of outputs of a set of parameter permutations from a dataframe `df`.
"""
function v_of_S(df::DataFrame, S::Union{Vector{Symbol}, Vector{String}}, y::Symbol)
    if isempty(S)
        return 0.0  # Var(E[Y|∅]) = 0 (E[Y|∅] is constant = overall mean)
    end
    g = combine(groupby(df, S), y => mean => :μ, nrow => :w)
    μ̄ = mean(df[!, y])
    # weighted population variance of group means around overall mean
    return sum(g.w .* (g.μ .- μ̄).^2) / sum(g.w)
end

"""
    shapley_main_and_interactions(df::DataFrame, variables::Union{Vector{Symbol}, Vector{String}}, output_col::Symbol)

Calculate Shapley values over sets of variables given permutations in `df` with output column
`y`. This function calculates main effect and pairwise interactions of variables.
"""
function shapley_main_and_interactions(df::DataFrame, variables::Union{Vector{Symbol}, Vector{String}}, output_col::Symbol)
    N = length(variables)
    total_var = Statistics.var(df[:, output_col])

    main_effects = Dict{Symbol,Float64}()
    interactions = Dict{Tuple{Symbol,Symbol},Float64}()

    # --- Main effects ---
    for var in variables
        φ_i = 0.0
        others = setdiff(variables, [var])
        for k in 0:length(others)
            for S in combinations(others, k)
                k = length(S)
                weight = factorial(k) * factorial(N-k-1) / factorial(N)
                v_S =  v_of_S(df, S, output_col)
                v_Si =  v_of_S(df, [S...; var], output_col)
                φ_i += weight * (v_Si - v_S)
            end
        end
        main_effects[var] = φ_i / total_var
    end

    # --- Pairwise interactions ---
    for var_i in variables
        j_vars = setdiff(variables, [var_i])
        for var_j in j_vars
            φ_ij = 0.0
            others = setdiff(variables, [var_i, var_j])
            for k in 0:length(others)
                for S in combinations(others, k)
                    k = length(S)
                    weight = factorial(k) * factorial(N-k-2) / factorial(N)
                    v_S   =  v_of_S(df, S, output_col)
                    v_Si  =  v_of_S(df, [S...; var_i], output_col)
                    v_Sj  =  v_of_S(df, [S...; var_j], output_col)
                    v_Sij =  v_of_S(df, [S...; var_i; var_j], output_col)
                    φ_ij += weight * (v_Sij - v_Si - v_Sj + v_S)
                end
            end
            interactions[(var_i, var_j)] = φ_ij / total_var
        end
    end

    return (main_effects = main_effects, interactions = interactions)
end

shap_effects = [
    shapley_main_and_interactions(df, Symbol.(input_variables), output) 
    for output in outputs]
main_effects = Dict(zip(
    outputs, 
    [effects.main_effects for effects in shap_effects]
))
variables_x_outputs = vcat(collect(Iterators.product(
    outputs,
    Symbol.(input_variables)
))...)

shap_effects_long = DataFrame(variable = last.(variables_x_outputs), output = first.(variables_x_outputs))
shap_effects_long.main_effect .= 0.0
for row in eachrow(shap_effects_long)
    row.main_effect = main_effects[row.output][row.variable]
end

output_plot = data(shap_effects_long) * mapping(:variable, :main_effect, row=:output) * visual(BarPlot, direction=:x)
fig = draw(output_plot)
fig



# Shapley effects for decade permutations

within_esm_ssp_perms = CSV.read("../outputs/within_esm_ssp_decade_permutations.csv", DataFrame)
input_variables = names(within_esm_ssp_perms)[(names(within_esm_ssp_perms) .!= "perm_id") .& (names(within_esm_ssp_perms) .!= "Column1")]
decade_encoding = Dict("2010-2019" => 1, "2030-2039" => 2, "2060-2069" => 3)

# For 2010-2019:
esm_ssp_perms_CNRM_SSP370 = within_esm_ssp_perms[contains.(within_esm_ssp_perms.light, "CNRM-ssp370"), :]
for var in input_variables
    esm_ssp_perms_CNRM_SSP370[!, var] = [find_encoded_val(x, decade_encoding, r"(-[[:upper:]]{4}-[[:lower:]]{3}\d{3})") for x in esm_ssp_perms_CNRM_SSP370[:, var]]
end

result_files = 1:200
result_files = ["../outputs/across_decade_permutations/model_outputs_perm_$(x).parq" for x in result_files]
df = esm_ssp_perms_CNRM_SSP370[1:200, :]
df.annual_surface_phyt = [extract_annual_results(filename, "Description", "Surface_layer_phytoplankton", "Model_annual_mean") for filename in result_files]
df.annual_dem_fish = [extract_annual_results(filename, "Description", "Demersal_fish", "Model_annual_mean") for filename in result_files]
df.annual_plank_fish = [extract_annual_results(filename, "Description", "Planktivorous_fish", "Model_annual_mean") for filename in result_files]
df.annual_nit_mass = [extract_annual_results(filename, "Description", "Total_nitrogen_mass", "Model_annual_mean") for filename in result_files]
df.annual_birds = [extract_annual_results(filename, "Description", "Birds", "Model_annual_mean") for filename in result_files]
outputs = [
    :annual_surface_phyt, 
    :annual_dem_fish, 
    :annual_plank_fish, 
    :annual_nit_mass,
    :annual_birds
]
shap_effects = [
    shapley_main_and_interactions(df, Symbol.(input_variables), output) 
    for output in outputs]
main_effects = Dict(zip(
    outputs, 
    [effects.main_effects for effects in shap_effects]
))
variables_x_outputs = vcat(collect(Iterators.product(
    outputs,
    Symbol.(input_variables)
))...)

shap_effects_long = DataFrame(variable = last.(variables_x_outputs), output = first.(variables_x_outputs))
shap_effects_long.main_effect .= 0.0
for row in eachrow(shap_effects_long)
    row.main_effect = main_effects[row.output][row.variable]
end

output_plot = data(shap_effects_long) * mapping(:variable, :main_effect, row=:output) * visual(BarPlot, direction=:x)
fig = draw(output_plot)
fig


# PCA analysis
result_dfs = [DataFrame(read_parquet(filename)) for filename in result_files] 
for (r, result_df) in enumerate(result_dfs)
    result_df[!, :perm_id] .= r
end
result_df = vcat(result_dfs...)
result_df = unstack(result_df, :perm_id, :Description, :Model_annual_mean)

result_df = result_df[:, (std.(eachcol(result_df)) .!= 0.0)]

annual_outputs = names(result_df)[names(result_df) .!= "perm_id"]
for col in annual_outputs
    μ, σ = rescale!(result_df[!, col]; obsdim=1)
end
# Remove any columns that only contain a constant value

heatmap(Matrix(result_df[:, annual_outputs]))

M_pca = fit(PCA, Matrix{Float64}(Matrix(result_df[:, annual_outputs])'); maxoutdim=2)
y_pca = predict(M_pca, Matrix{Float64}(Matrix(result_df[:, annual_outputs])'))
pca_val_df = DataFrame(perm_id = result_df.perm_id, PC1 = y_pca[1, :], PC2 = y_pca[2, :])
draw(data(pca_val_df) * mapping(:PC1, :PC2, color=:perm_id) * visual(Scatter))

M_kpca = fit(KernelPCA, Matrix{Float64}(Matrix(result_df[:, annual_outputs])'); maxoutdim=2)
y_kpca = predict(M_kpca, Matrix{Float64}(Matrix(result_df[:, annual_outputs])'))
kpca_val_df = DataFrame(perm_id = result_df.perm_id, PC1 = y_kpca[1, :], PC2 = y_kpca[2, :])
kpca_val_df.decade .= ""
for row in eachrow(kpca_val_df)
    perm = row.perm_id
    decades = within_esm_ssp_perms[within_esm_ssp_perms.perm_id .== perm, input_variables]
    decades = [first(decades[:, var]) for var in input_variables]
    decades = [replace(entry, r"(-[[:upper:]]{4}-[[:lower:]]{3}\d{3})" => "") for entry in decades]
    row.decade = mode(decades)
end
draw(data(kpca_val_df) * mapping(:PC1, :PC2, color=:decade) * visual(Scatter))

