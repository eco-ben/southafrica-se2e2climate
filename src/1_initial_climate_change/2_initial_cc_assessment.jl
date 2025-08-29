using CSV, DataFrames
using Parquet
using DataFrames
using Statistics
using MultivariateStats
using AlgebraOfGraphics
using GLMakie

master = CSV.read("../outputs/master_forcings_South_Africa_MA.csv", DataFrame)

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
result_df_lines = result_df[result_df.Description .∈ [output_variables], :]
result_df_lines = result_df_lines[result_df_lines.Description .∈ [result_df_lines.Description[contains.(result_df_lines.Description, ["fish"])]], :]
result_df_lines.decade = [first(match(r"(\d{4}-\d{4})", var).captures) for var in result_df_lines.variant]
result_df_lines.ESM_SSP = [first(match(r"([[:upper:]]{4}-[[:lower:]]{3}\d{3})", var).captures) for var in result_df_lines.variant]
biomass_timeseries = data(result_df_lines) * mapping(:decade, :Model_annual_mean, color=:ESM_SSP, layout=:Description) * visual(Lines)
biomass_ts_fig = draw(biomass_timeseries; facet=(; linkyaxes=:none), figure=(;fontsize=8, size = (900, 600)), legend=(;position=:bottom))
save("../figs/initial_cc_assessment/biomass_timeseries.png", biomass_ts_fig)

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
pca_line_fig = draw(scat + line; figure=(;fontsize=8, size=(600, 500)))
save("../figs/initial_cc_assessment/biomass_pca.png", pca_line_fig)

L = projection(M_pca)
ev = principalvars(M_pca) ./ sum(principalvars(M_pca))  # variance explained
corr_circle = L .* sqrt.(ev[1:2]')                      # scale loadings

scaling_factor = maximum(abs, y_pca) / 1.5   # heuristic for visibility
L_scaled = corr_circle .* scaling_factor
corr_circle = L_scaled

fig = pca_line_fig.figure
# ax2 = Axis(fig[2, 1], xlabel="PC1", ylabel="PC2")

# draw unit circle (correlation circle)
# θ = range(0, 2π, length=200)
# lines!(ax, cos.(θ), sin.(θ), color=:gray, linestyle=:dash)

# Draw arrows for variables
for i in 1:Base.size(corr_circle, 1)
    arrows2d!(fig.layout[1,1], [mean(y_pca[1, :])], [mean(y_pca[2, :])], [corr_circle[i,1]], [corr_circle[i,2]], 
            shaftwidth=2, color=:blue, alpha=0.6)
    text!(fig.layout[1,1], corr_circle[i,1], corr_circle[i,2], text=annual_outputs[i], align=(:left, :bottom), alpha=0.5)
end

fig
