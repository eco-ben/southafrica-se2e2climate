using CSV, DataFrames, StatsBase, StatsPlots
using Statistics

physics = CSV.read("../../../StrathE2E_workspace/Models/South_Africa_MA/2010-2019-CNRM-ssp126/Driving/physics_SOUTH_AFRICA_MA_2010-2019-CNRM-ssp126.csv", DataFrame)
chemistry = CSV.read("../../../StrathE2E_workspace/Models/South_Africa_MA/2010-2019-CNRM-ssp126/Driving/chemistry_SOUTH_AFRICA_MA_2010-2019-CNRM-ssp126.csv", DataFrame)
forcing = hcat(physics, chemistry[:, Not(:Month)])
selected_cols = names(forcing)[names(forcing) .∉ [[
    "Month"; 
    names(forcing)[contains.(names(forcing), ["DO"])]; 
    names(forcing)[contains.(names(forcing), ["pdist"])];
    names(forcing)[contains.(names(forcing), ["SPM"])]
]]]
forcing = forcing[:, selected_cols]

@df forcing_10_19_CNRM_126 cornerplot(cols(8:14))

cor_mat = corspearman(Matrix(forcing))

fig = Figure()
ax = Axis(
    fig[1,1],
    xticks = (1:ncol(forcing), names(forcing)),
    yticks = (1:ncol(forcing), names(forcing)),
    xticklabelrotation = 45
)
hm = GLMakie.heatmap!(cor_mat, colormap=:balance)
GLMakie.Colorbar(fig[1,2], hm)

variable_groups = Dict(
    "temp" => names(forcing)[contains.(names(forcing), ["temp"])],
    "light" => "SLight",
    "nutrient_concentrations" => names(forcing)[
        (contains.(names(forcing), [r"nitrate$"]) .| 
        contains.(names(forcing), [r"ammonia$"]) .| 
        contains.(names(forcing), [r"phyt$"]) .| 
        contains.(names(forcing), [r"detritus$"])) .&
        .!contains.(names(forcing), [r"^RIV"])
    ],
    "river_inputs" => names(forcing)[contains.(lowercase.(names(forcing)), [r"^riv"])],
    "nutrient_flux" => names(forcing)[contains.(names(forcing), [r"flux$"])],
    "boundary_flow" => names(forcing)[
        contains.(names(forcing), [r"OUT$"]) .|
        contains.(names(forcing), [r"IN$"])
    ],
    "vertical_flow" => ["log10Kvert", "mixLscale", "D_SO_upwelling", "SO_D_downwelling"],
)