## This will activate the environment in the code sub-folder
import Pkg
Pkg.activate("code")

# Generally, it is considered good practice to work in an environment, as it
# makes the management of dependencies a lot easier, and the project becomes
# more reproducible

## We then load the packages
using LinearAlgebra, Statistics
using EcologicalNetworks
using StatsBase
using Gadfly
import Cairo, Fontconfig
using Colors
using DataFrames, Query
using GLM, ANOVA

## Import the functions and methods we need
include(joinpath(pwd(), "code", "lib", "main.jl"))
include(joinpath(pwd(), "code", "lib", "extinctions.jl"))

## Using the Web Of Life dataset
networks = DataFrame(
    ID = String[],
    type = Symbol[],
    latitude = Union{Missing,Float64}[],
    longitude = Union{Missing,Float64}[],
    entropy = Float64[],
    defficiency = Float64[],
    richness = Int64[],
    s1 = Int64[],
    s2 = Int64[],
    nestedness = Float64[],
    spectral_radius = Float64[],
    connectance = Float64[]
    )

for wol in web_of_life()
    N = simplify(convert(BipartiteNetwork, web_of_life(wol.ID)))
    if richness(N) <= 300
        D = Dict{Symbol, Any}()
        D[:ID] = wol.ID
        D[:latitude] = wol.Latitude == "" ? missing : wol.Latitude
        D[:longitude] = wol.Longitude == "" ? missing : wol.Longitude
        D[:type] = Symbol(wol.Type_of_interactions)
        D[:entropy] = svd_entropy(N)
        D[:nestedness] = η(N)
        D[:spectral_radius] = ρ(N)
        D[:richness] = richness(N)
        D[:connectance] = connectance(N)
        D[:defficiency] = ((maxrank(N) - rank(N)) / maxrank(N))
        D[:s1] = richness(N, dims=1)
        D[:s2] = richness(N, dims=2)
        push!(networks, D)
    end
end

mmin = (x) -> minimum(filter(!ismissing, x))
mmax = (x) -> maximum(filter(!ismissing, x))

gdf = groupby(networks, :type)
summ = combine(gdf,
    :ID => length => :n,
    :latitude => mmin => :minlat,
    :latitude => mmax => :maxlat,
    :s1 => mean => :rtop,
    :s2 => mean => :rdown
)
sort!(summ, :n, rev=true)

# Colour palette
colour_palette = Scale.color_discrete_manual(
    colorant"#648FFF",
    colorant"#785EF0",
    colorant"#DC267F",
    colorant"#FE6100",
    colorant"#FFB000"
    )

## Theme for the plots

paper_theme = Theme(
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
Gadfly.set_default_plot_size(18cm, 12cm)

## ANOVA for interaction types

#remove Networks types with low n
q1 = @from i in networks begin
            @where i.type != Symbol("Plant-Ant")
            @where i.type != Symbol("Plant-Herbivore")
            @select {i.entropy, i.type}
            @collect DataFrame
       end

model = fit(LinearModel,
            @formula(entropy ~  type),
            q1,
            contrasts = Dict(:type => EffectsCoding()))
anova(model)

#=TODO
add Tukey HSD=#

## Interaction Type vs Rank & Entropy
draw(
    PNG(joinpath("figures", "interactiontype_v_entropy.png"), dpi=300),
    plot(networks,
        x=:type, y=:entropy,
        color=:type,
        Geom.beeswarm,
        Guide.xlabel("Interaction Type"), Guide.ylabel("Entropy"),
        Scale.x_discrete(order=[1,5,2,3,4])
    )
)

## Entropy and relative rank defficiency
draw(
    PNG(joinpath("figures", "entropy_v_rank.png"), dpi = 300),
    plot(networks,
        x=:defficiency, y=:entropy,
        color=:type,
        Geom.point,
        Guide.xlabel("Relative rank defficiency"), Guide.ylabel("Entropy")
    )
)

## Size vs Rank & Entropy

#= TODO
* Change labels on Y axes (this could come down to changing things in the DataFrame?)
* Labels on subplots e.g. A/B (see below for an attempt) using layers =#

ms = Dict(:defficiency=>"Relative rank defficiency", :entropy=>"SVD entropy")

draw(
    PNG(joinpath("figures", "size_v_rankentropy.png"), dpi=300),
    plot(
        stack(networks, [:defficiency, :entropy], variable_name =:measure, value_name=:value),
        x=:richness, y=:value,
        color=:type, ygroup =:measure,
        #Scale.ygroup(labels=i->f(i), levels=collect(keys(ms))), # this should work, it doesn't
        Geom.subplot_grid(
            Geom.point, free_y_axis=true
        ),
        Guide.xlabel("Richness"), Guide.ylabel(nothing)
    )
)

#= 1st pass at labelling sub plots... its not going well...
plot(stack(Outputs, [:RankDefficiencyRel, :Entropy], variable_name =:measure, value_name=:value),
ygroup =:measure,
Geom.subplot_grid(
layer(
x=:Richness, y=:value, color=:InteractionType, Geom.point),
layer(DataFrame(x = [0.1,0.1],
y = [0.51,1.1],
label = ["A", "B"],
group = ["RankDefficiencyRel", "Entropy"]),
x=:x, y=:y, label =:label, ygroup =:group, Geom.label),
free_y_axis=true),
alpha = [0.6], Guide.xlabel("Richness"), Guide.ylabel(nothing))=#

## Other measures of networks vs Rank & Entropy

#= TODO
* Change labels on X axes (this could come down to changing things in the DataFrame?)
* Labels on subplots e.g. A/B (see below for an attempt) using layers =#

draw(
    PNG(joinpath("figures", "others_v_entropy.png"), dpi = 300),
    plot(
        stack(
            stack(networks, [:entropy], [:nestedness, :spectral_radius, :connectance, :type], variable_name =:measure, value_name=:value),
            [:nestedness, :spectral_radius, :connectance],
        variable_name =:method, value_name=:metric),
        x=:metric, y=:value,
        color=:type, ygroup =:measure, xgroup =:method, Geom.subplot_grid(Geom.point, free_x_axis=true),
        alpha = [0.8], Guide.xlabel(nothing), Guide.ylabel(nothing)
    )
)

## Resilience vs Rank & Entropy
f_random = _order_species_for_removal();
f_increasing = _order_species_for_removal(true);
f_decreasing = _order_species_for_removal(false);

AUC = DataFrame(
    ID = String[],
    type = Symbol[],
    entropy = Float64[],
    Random_all = Float64[],
    Increasing_all = Float64[],
    Decreasing_all = Float64[],
    Random_1 = Float64[],
    Increasing_1 = Float64[],
    Decreasing_1 = Float64[],
    Random_2 = Float64[],
    Increasing_2 = Float64[],
    Decreasing_2 = Float64[]
    )

for wol in web_of_life()
    N = simplify(convert(BipartiteNetwork, web_of_life(wol.ID)))
    if richness(N) <= 200
        D = Dict{Symbol, Any}()
        D[:ID] = wol.ID
        D[:type] = Symbol(wol.Type_of_interactions)
        D[:entropy] = svd_entropy(N)
        D[:Random_all] = extinction_robustness(extinction(N, f_random))
        D[:Increasing_all] = extinction_robustness(extinction(N, f_increasing))
        D[:Decreasing_all] = extinction_robustness(extinction(N, f_decreasing))
        D[:Random_1] = extinction_robustness(extinction(N, f_random, dims = 1), dims = 2)
        D[:Increasing_1] = extinction_robustness(extinction(N, f_increasing, dims = 1), dims = 2)
        D[:Decreasing_1] = extinction_robustness(extinction(N, f_decreasing, dims = 1), dims = 2)
        D[:Random_2] = extinction_robustness(extinction(N, f_random, dims = 2), dims = 1)
        D[:Increasing_2] = extinction_robustness(extinction(N, f_increasing, dims = 2), dims = 1)
        D[:Decreasing_2] = extinction_robustness(extinction(N, f_decreasing, dims = 2), dims = 1)
        push!(AUC, D)
    end
end

#Melt data into long format
AUC = stack(AUC, [:Random_all, :Random_1, :Random_2, :Decreasing_all, :Decreasing_1,
:Decreasing_2, :Increasing_all, :Increasing_1, :Increasing_2],
variable_name =:measure, value_name=:value);
#Split out by dimension and extinction mechanism
AUC = hcat(AUC, DataFrame([(Type=a,Dimension=b) for (a,b) in split.(AUC.measure, "_")]))

#= TODO
* Change labels on X axes (this could come down to changing things in the DataFrame?)
* Labels on subplots e.g. A/B (see below for an attempt) using layers
* decide how we want to refer to 'mechanism' and 'dimension' =#
draw(
    PNG(joinpath("figures", "entropy_v_AUCall.png"), dpi = 300),
    plot(AUC,
        x =:value, y =:Entropy,
        color=:InteractionType, ygroup =:Type, xgroup =:Dimension,
        Geom.subplot_grid(Geom.point, free_y_axis=true),
        Guide.xlabel("Resilience"), Guide.ylabel("Extinction mechanism")
    )
)

## End of script
