import Pkg
Pkg.activate("code")

using LinearAlgebra, Statistics
using EcologicalNetworks
using StatsBase
using Gadfly
import Cairo, Fontconfig
using Colors
using DataFrames

colour_palette = Scale.color_discrete_manual(
    colorant"#648FFF",
    colorant"#785EF0",
    colorant"#DC267F",
    colorant"#FE6100",
    colorant"#FFB000"
    )

## Theme for the plots

paper_theme = Gadfly.Theme(
    panel_stroke = colorant"#000000",
    panel_fill = colorant"#ffffff",
    panel_line_width = 0.5mm,
    plot_padding = [5mm],
    grid_color = colorant"#dcdedf",
    grid_line_style = :solid,
    major_label_color = colorant"#222222",
    minor_label_color = colorant"#000000",
    key_label_color = colorant"#000000",
    key_title_color = colorant"#000000",
    guide_title_position = :center,
    colorkey_swatch_shape = :circle,
    key_position = :top,
    key_title_font_size = 0mm,
    discrete_color_scale = colour_palette
)
Gadfly.push_theme(paper_theme)
Gadfly.set_default_plot_size(16cm, 16cm)

include(joinpath(pwd(), "code", "lib", "main.jl"))
include(joinpath(pwd(), "code", "lib", "extinctions.jl"))
include(joinpath(pwd(), "code", "lib", "annealing.jl"))

function generate_sample(n, m, r)
    eq = isequal(size(n))
    return filter(k -> eq(size(k)), simplify.(rand(m(n), r)))
end

# Analysis
res = DataFrame(
    id = String[],
    type = String[],
    entropy = Float64[],
    richness = Int64[],
    links = Int64[],
    connectance = Float64[],
    model = [],
    zscore = Float64[]
)

n_small = filter(n -> n.Species <= 200, web_of_life())
filter!(n -> n.Type_of_interactions ∈ ["Host-Parasite", "Pollination", "Seed Dispersal"], n_small)
for n in n_small
    N = convert(BipartiteNetwork, web_of_life(n.ID))
    for m in [null1, null2]
        ℕ = svd_entropy.(generate_sample(N, m, 999))
        z = (svd_entropy(N) - mean(ℕ))/std(ℕ)
        push!(res,
            (
                n.ID,
                n.Type_of_interactions,
                svd_entropy(N),
                EcologicalNetworks.richness(N),
                links(N),
                connectance(N),
                m,
                z
            )
        )
    end
end

f = (x) -> 1.0/(1.0+exp(-x))
successes = res[findall(.!isnan.(res.zscore)),:]
successes = successes[findall(successes.zscore .>= -20.0),:]
successes.score = f.(successes.zscore)

replace!(successes.model, null1 => "Type 1");
replace!(successes.model, null2 => "Type 2");

p1 = Gadfly.plot(successes,
    x = :zscore,
    color = :type,
    ygroup = :model,
    xgroup = :type,
    Geom.subplot_grid(
        Geom.histogram,
        Guide.xticks(orientation=:horizontal)
    ),
)

p2 = Gadfly.plot(successes,
    x = :richness,
    y = :score,
    color = :type,
    ygroup = :model,
    xgroup = :type,
    Geom.subplot_grid(
        Geom.point,
        Scale.x_log2,
        Guide.xticks(orientation=:horizontal)
    ),
)

draw(
    PNG(joinpath("figures", "nullmodel_histogram.png"), dpi = 300),
    p1
)

draw(
    PNG(joinpath("figures", "nullmodel_richness.png"), dpi = 300),
    p2
)
