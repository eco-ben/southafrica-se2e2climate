using CSV, DataFrames
using Parquet
using DataFrames
using Statistics
using MultivariateStats
using AlgebraOfGraphics
using GLMakie

master = CSV.read("../outputs/master_forcings_South_Africa_MA.csv", DataFrame)

function find_encoded_val(full_code_value, encoding_scheme, disregarded_pattern)
    split_code_value = replace(full_code_value, disregarded_pattern => "")

    return encoding_scheme[split_code_value]
end

function extract_annual_results(filename, id_column, id, value_column)
    df = DataFrame(read_parquet(filename))

    return first(df[df[:, id_column] .== id, value_column])
end

result_files = readdir("../outputs/initial_runs/", join=true)
result_files = result_files[contains.(result_files, ["final_biomass"])]
variants = [first(match(r"(\d{4}-\d{4}-[[:upper:]]{4}-[[:lower:]]{3}\d{3})", fn).captures) for fn in result_files]

result_dfs = [DataFrame(read_parquet(filename)) for filename in result_files] 
for (r, result_df) in enumerate(result_dfs)
    result_df[!, :variant] .= variants[r]
end
result_df = vcat(result_dfs...)

output_variables = unique([
    "variant";
    result_df.Description[contains.(result_df.Description, ["detritus"])];
    result_df.Description[contains.(result_df.Description, ["plankton"])];
    result_df.Description[contains.(result_df.Description, ["fish"])];
    result_df.Description[contains.(result_df.Description, ["Benthos"])];
    ["Birds", "Pinnipeds", "Cetaceans"]
])

# Plot timeseries of interesting biomass changes under climate change
result_df_lines = result_df_lines[result_df_lines.Description .∈ [output_variables], :]
result_df_lines = result_df_lines[result_df_lines.Description .∈ [result_df_lines.Description[contains.(result_df_lines.Description, ["fish"])]], :]
result_df_lines.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in result_df_lines.variant]
result_df_lines.ESM_SSP = [first(match(r"([[:upper:]]{4}-[[:lower:]]{3}\d{3})", var).captures) for var in result_df_lines.variant]
biomass_timeseries = data(result_df_lines) * mapping(:decade, :Model_annual_mean, color=:ESM_SSP, layout=:Description) * visual(Lines)
draw(biomass_timeseries; facet=(; linkyaxes=:none), figure=(;fontsize=9))

# PCA modelling of the biomass results
result_df_wide = unstack(result_df, :variant, :Description, :Model_annual_mean)
result_df_wide = result_df_wide[:, output_variables]

annual_outputs = required_cols[required_cols .!= "variant"]
for col in annual_outputs
    μ, σ = rescale!(result_df_wide[!, col]; obsdim=1)
end
# Remove any columns that only contain a constant value

heatmap(Matrix(result_df_wide[:, annual_outputs]))

M_pca = fit(PCA, Matrix{Float64}(Matrix(result_df_wide[:, annual_outputs])'); maxoutdim=2)
y_pca = predict(M_pca, Matrix{Float64}(Matrix(result_df_wide[:, annual_outputs])'))
pca_val_df = DataFrame(variant = result_df_wide.variant, PC1 = y_pca[1, :], PC2 = y_pca[2, :])

pca_val_df.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in pca_val_df.variant]
pca_val_df.ESM = [first(match(r"([[:upper:]]{4})", var).captures) for var in pca_val_df.variant]
pca_val_df.SSP = [first(match(r"([[:lower:]]{3}\d{3})", var).captures) for var in pca_val_df.variant]
pca_val_df.ESM_SSP = [first(match(r"([[:upper:]]{4}-[[:lower:]]{3}\d{3})", var).captures) for var in pca_val_df.variant]

scat = data(pca_val_df) * mapping(:PC1, :PC2, color=:ESM_SSP, marker=:decade) * visual(Scatter)
line = data(sort(pca_val_df, :decade)) * mapping(:PC1, :PC2, color=:ESM_SSP) * visual(Lines)
draw(scat + line)