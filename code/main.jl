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
                    SpectralRadius = ρ.(Bs),
                    InteractionType = [y[:Type_of_interactions] for y in web_of_life()]) begin
            @where i.Richness < 200
            @select {i.Entropy, i.RankDefficiencyRel, i.Richness, i.Nestedness, i.SpectralRadius, i.InteractionType}
            @collect DataFrame
       end

## Interaction Type vs Rank & Entropy

draw(PNG(joinpath("figures", "interactiontype_v_entropy.png"), dpi=300),
    plot(Outputs,
        x=:InteractionType, y=:Entropy,
        color=:InteractionType,
        Geom.beeswarm, alpha = [0.3],
        Guide.xlabel("Interaction Type"), Guide.ylabel("Entropy")))

## Here we could plot Entropy and various measures of Rank

draw(PNG(joinpath("figures", "entropy_v_rank.png"), 20cm, 20cm, dpi = 300),
plot(Outputs,
    x=:RankDefficiencyRel, y=:Entropy,
    color=:InteractionType, alpha = [0.6], Guide.xlabel("Relative rank defficiency"), Guide.ylabel("Entropy")))

## Size vs Rank & Entropy

draw(PNG(joinpath("figures", "size_v_rank&entropy.png"), 15cm, 20cm, dpi=300),
    plot(stack(Outputs, [:RankDefficiencyRel, :Entropy], variable_name =:measure, value_name=:value),
        x=:Richness, y=:value,
        color=:InteractionType, ygroup =:measure, Geom.subplot_grid(Geom.point, free_y_axis=true),
        alpha = [0.6], Guide.xlabel("Richness"), Guide.ylabel(nothing)))

## Other measures of networks vs Rank & Entropy

draw(PNG(joinpath("figures", "others_v_entropy.png"), dpi = 300),
    plot(stack(stack(Outputs, [:Entropy], [:Nestedness, :SpectralRadius, :InteractionType],
        variable_name =:measure, value_name=:value), [:Nestedness, :SpectralRadius],
        variable_name =:method, value_name=:metric),
        x=:metric, y=:value,
        color=:InteractionType, ygroup =:measure, xgroup =:method, Geom.subplot_grid(Geom.point, free_y_axis=true),
        alpha = [0.6], Guide.xlabel(nothing), Guide.ylabel(nothing)))

## Resilience vs Rank & Entropy

f_random = _order_species_for_removal();
f_increasing = _order_species_for_removal(true);
f_decreasing = _order_species_for_removal(false);

AUC = DataFrame(Random_all = [extinction_robustness(extinction(B, f_random)) for B in Bs],
                Increasing_all = [extinction_robustness(extinction(B, f_increasing)) for B in Bs],
                Decreasing_all = [extinction_robustness(extinction(B, f_decreasing)) for B in Bs],
                Random_1 = [extinction_robustness(extinction(B, f_random, dims = 1), dims = 2) for B in Bs],
                Increasing_1 = [extinction_robustness(extinction(B, f_increasing, dims = 1), dims = 2) for B in Bs],
                Decreasing_1 = [extinction_robustness(extinction(B, f_decreasing, dims = 1), dims = 2) for B in Bs],
                Random_2 = [extinction_robustness(extinction(B, f_random, dims = 2), dims = 1) for B in Bs],
                Increasing_2 = [extinction_robustness(extinction(B, f_increasing, dims = 2), dims = 1) for B in Bs],
                Decreasing_2 = [extinction_robustness(extinction(B, f_decreasing, dims = 2), dims = 1) for B in Bs],
                Entropy = svd_entropy.(Bs),
                InteractionType = [y[:Type_of_interactions] for y in web_of_life()])
AUC = stack(AUC, [:Random_all, :Random_1, :Random_2, :Decreasing_all, :Decreasing_1,
            :Decreasing_2, :Increasing_all, :Increasing_1, :Increasing_2],
            variable_name =:measure, value_name=:value);
AUC = hcat(AUC, DataFrame([(Type=a,Dimension=b) for (a,b) in split.(AUC.measure, "_")]))

draw(PNG(joinpath("figures", "entropy_v_AUCall.png"), dpi = 300),
        plot(AUC,
            x =:value, y =:Entropy,
            color=:InteractionType, ygroup =:Type, xgroup =:Dimension,
            Geom.subplot_grid(Geom.point, free_y_axis=true),
            alpha = [0.6], Guide.xlabel("Dimension"), Guide.ylabel("Extinction method")))

## End of script
