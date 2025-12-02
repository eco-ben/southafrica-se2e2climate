using Rasters, GLMakie
using ColorSchemes
using Parquet
using DataFrames
using Dates
import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import NCDatasets

inch = 96

land = GDF.read("../../../Spatial Data/ne_10m_land/ne_10m_land.shp")[1, :]
habitats = GDF.read("../data/Habitats.gpkg")
habitats.class = habitats.Shore .* " - " .* habitats.Habitat

class_colors = Dict(
    "Inshore - rock" => "#40333C",
    "Inshore - mud" => "#284481",
    "Inshore - sand" => "#9097CC",
    "Inshore - gravel" => "#4A8FA1",
    "Offshore - rock" => "#d7c288",
    "Offshore - mud" => "#ffb700",
    "Offshore - sand" => "#FFD25F",
    "Offshore - gravel" => "#ffedbd"
)
habitats = habitats[indexin(collect(keys(class_colors)), habitats.class), :]

fig = Figure()
ax = Axis(
    fig[1,1],
    limits = ((14, 23), (-37.5, -28.5)),
    backgroundcolor=(:lightblue, 0.3)
)
poly!(ax, land.geometry, color = (:gray, 0.8))
poly!(ax, habitats.geometry, color = 1:8, colormap=collect(values(class_colors)))

legend_entries = [PolyElement(; color = class_colors[class]) for class in habitats.class]
Legend(fig[2,1], legend_entries, habitats.class, orientation=:horizontal, nbanks = 2)
save("../figs/habitat_map.png", fig, px_per_unit = 300/inch)


# Plot map with heatmap of indicative fishing
gfw = DataFrame(read_parquet("../../../Spatial Data/fishing_effort_data/Global_fishing_watch/fleet-daily.parq"))
gfw.date = Date.(gfw.date)
gfw = gfw[
    (gfw.cell_ll_lat .< -28) .&
    (gfw.cell_ll_lat .> -38) .&
    (gfw.cell_ll_lon .< 30) .&
    (gfw.cell_ll_lon .> 14) .&
    (gfw.date .< Date("2020-01-01")),
:]

gfw_sum = DataFrames.combine(groupby(gfw, [:cell_ll_lat, :cell_ll_lon]), :hours => sum)
scatter(gfw_sum.cell_ll_lon, gfw_sum.cell_ll_lat, color = gfw_sum.hours_sum, alpha = 0.2)
lon_range = extrema(gfw_sum.cell_ll_lon)
lat_range = extrema(gfw_sum.cell_ll_lat)

lon, lat = X(first(lon_range):0.01:last(lon_range)), Y(first(lat_range):0.01:last(lat_range))
gfw_ras = Raster(Array{Union{Missing, Float64}}(undef, lon, lat), dims = (lon, lat))

for (lon, lat) in unique(zip(gfw_sum.cell_ll_lon, gfw_sum.cell_ll_lat))
    gfw_ras[X = At(lon), Y = At(lat)] = first(gfw_sum[(gfw_sum.cell_ll_lon .== lon) .& (gfw_sum.cell_ll_lat .== lat), :].hours_sum)
end

gfw_ras = ifelse.((.!ismissing.(gfw_ras)) .&(gfw_ras .< 10.0), missing, gfw_ras) # Remove cells that have less than 10h over the years
gfw_ras_log = log10.(gfw_ras)

fig = Figure()
ax = Axis(
    fig[1,1],
    limits = ((13, 30), (-37.5, -27.5))
)
poly!(ax, land.geometry, color = (:gray, 0.8))
hm = heatmap!(ax, gfw_ras_log)
poly!(ax, habitats.geometry, color = (:red, 0.2))

Colorbar(fig[1,2], hm, label = "Log10 total fishing hours (>10h)")
save("../figs/gfw_fishing_domain_map.png", fig, px_per_unit = 300/inch)


# Plot bathymetry and distance to shore map
sabounds = X(12 .. 40), Y(-39 .. -25) # Roughly cut out SA

dist = Raster("../../Shared Data/distance-from-shore.tif"; lazy=true)
dist = dist[sabounds...] |> trim
dist = ifelse.(((dist .== 0) .| (dist .> 20)), missing, dist)

bathy = Raster("../../Shared Data/GEBCO_2020.nc"; lazy=true)
bathy = bathy[sabounds...] |> trim
bathy = ifelse.(((bathy .>= 0) .| (bathy .< -800)), missing, bathy) # Remove cells that have less than 10h over the years
bathy = .-bathy

fig = Figure()
ax1 = Axis(
    fig[1,1],
    limits = ((14, 23), (-37.5, -28.5))
)

poly!(ax1, land.geometry, color = :gray)
hm1 = heatmap!(ax1, dist)
poly!(ax1, habitats.geometry, color = (:red, 0.2))
Colorbar(fig[1,0], hm1, label="Distance to shore (maximum 20km)", flipaxis=false)

ax2 = Axis(
    fig[1,2],
    limits = ((14, 23), (-37.5, -28.5))
)

poly!(ax2, land.geometry, color = :gray)
hm2 = heatmap!(ax2, bathy)
poly!(ax2, habitats.geometry, color = (:red, 0.2))
Colorbar(fig[1,3], hm2, label="Bathymetry (maximum 800m)")

save("../figs/domain_separation_map.png", fig, px_per_unit = 300/inch)