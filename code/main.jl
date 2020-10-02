## This will activate the environment in the code sub-folder
import Pkg
Pkg.activate("code")

# Generally, it is considered good practice to work in an environment, as it
# makes the management of dependencies a lot easier, and the project becomes
# more reproducible

## We then load the packages
using LinearAlgebra
using EcologicalNetworks
using Statistics
using StatsBase
using Gadfly
using Cairo
using Fontconfig
using DataFrames
using Query

## Import the functions and methods we need
include(joinpath(pwd(), "code", "lib", "main.jl"))
include(joinpath(pwd(), "code", "lib", "extinctions.jl"))

## Using the Web Of Life dataset
raw_wol = [web_of_life(x.ID) for x in web_of_life()]
Bs = convert.(BipartiteNetwork, raw_wol)

## Wrangle data
Outputs = @from i in DataFrame(Richness = richness.(Bs),
                    Entropy = svd_entropy.(Bs),
                    RankDefficiencyRel = ((maxrank.(Bs) .- rank.(Bs)) ./ maxrank.(Bs)),
                    Nestedness = η.(Bs),
                    SpectralRadiance = ρ.(Bs),
                    InteractionType = [y[:Type_of_interactions] for y in web_of_life()]) begin
            @where i.Richness < 200
            @select {i.Entropy, i.RankDefficiencyRel, i.Richness, i.Nestedness, i.SpectralRadiance, i.InteractionType}
            @collect DataFrame
       end

## Interaction Type vs Rank & Entropy

draw(PNG(joinpath("figures", "interactiontype_v_rank&entropy.png"), 15cm, 20cm, dpi=300),
    plot(stack(Outputs, [:RankDefficiencyRel, :Entropy], variable_name =:measure, value_name=:value),
        x=:InteractionType, y=:value,
        color=:InteractionType, ygroup =:measure, Geom.subplot_grid(Geom.boxplot, free_y_axis=true),
        Guide.xlabel("Interaction Type"), Guide.ylabel(nothing)))

## Size vs Rank & Entropy

draw(PNG(joinpath("figures", "size_v_rank.png"), 15cm, 20cm, dpi=300),
    plot(stack(Outputs, [:RankDefficiencyRel, :Entropy], variable_name =:measure, value_name=:value),
        x=:Richness, y=:value,
        color=:InteractionType, ygroup =:measure, Geom.subplot_grid(Geom.point, free_y_axis=true),
        alpha = [0.6], Guide.xlabel("Richness"), Guide.ylabel(nothing),
        Guide.yticks("RankDefficiencyRel", "Entropy")))

## Here we could plot Entropy and various measures of Rank

#=
draw(PNG(joinpath("figures", "entropy_v_rank.png"), 20cm, 20cm, dpi = 300),
    vstack(hstack(plot(Outputs, x=:Richness, y=[:RankDefficiencyRel, :Entropy], Geom.point),
                plot(Outputs, x=:Rank, y=:Entropy, Geom.point)),
            hstack(plot(Outputs, x=:RankDefficiency, y=:Entropy, Geom.point, Guide.xlabel("Rank Defficiency")),
                plot(Outputs, x=:RankDefficiencyRel, y=:Entropy, Geom.point, Guide.xlabel("Rank Defficiency (relative)")))))
=#

## Other measures of networks vs Rank & Entropy

draw(PNG(joinpath("figures", "entropy_v_others.png"), 20cm, 20cm, dpi = 300),
    plot(stack(stack(Outputs, [:RankDefficiencyRel, :Entropy], [:Nestedness, :SpectralRadiance, :InteractionType],
        variable_name =:measure, value_name=:value), [:Nestedness, :SpectralRadiance],
        variable_name =:method, value_name=:metric),
        x=:metric, y=:value,
        color=:InteractionType, ygroup =:measure, xgroup =:method, Geom.subplot_grid(Geom.point, free_y_axis=true),
        alpha = [0.6], Guide.xlabel(nothing), Guide.ylabel(nothing)))

## Resilience vs Rank & Entropy

f_random = _order_species_for_removal();
f_increasing = _order_species_for_removal(true);
f_decreasing = _order_species_for_removal(false);

AUC = DataFrame(Random = [extinction_robustness(extinction(B, f_random)) for B in Bs],
                Increasing = [extinction_robustness(extinction(B, f_increasing, dims = 1)) for B in Bs],
                Decreasing = [extinction_robustness(extinction(B, f_decreasing, dims = 1)) for B in Bs],
                Rank = ((maxrank.(Bs) .- rank.(Bs)) ./ maxrank.(Bs)),
                Entropy = svd_entropy.(Bs),
                InteractionType = [y[:Type_of_interactions] for y in web_of_life()])
AUC = stack(AUC, [:Rank, :Entropy], variable_name =:measure, value_name=:value)
AUC = stack(AUC, [:Random, :Increasing, :Decreasing], variable_name =:method, value_name=:auc)

draw(PNG(joinpath("figures", "rank&entropy_v_AUC.png"), 30cm, 20cm, dpi = 300),
    plot(AUC, x =:auc, y =:value,
        color=:InteractionType, ygroup =:measure, xgroup =:method, Geom.subplot_grid(Geom.point, free_y_axis=true),
        Theme(key_position=:none), alpha = [0.6], Guide.xlabel(nothing), Guide.ylabel(nothing)))


## End of script
