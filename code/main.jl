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
using Colors

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
            #remove networks richness >200
            @where i.Richness < 200
            @select {i.Entropy, i.RankDefficiencyRel, i.Richness, i.Nestedness, i.SpectralRadius, i.InteractionType}
            @collect DataFrame
       end

#specify colour palette
ColourPalette = Scale.color_discrete_manual(colorant"#648FFF",
                            colorant"#785EF0",
                            colorant"#DC267F",
                            colorant"#FE6100",
                            colorant"#FFB000")

## Interaction Type vs Rank & Entropy

draw(PNG(joinpath("figures", "interactiontype_v_entropy.png"), dpi=300),
    plot(Outputs,
        x=:InteractionType, y=:Entropy,
        color=:InteractionType,
        Geom.beeswarm, alpha = [0.3],
        Guide.xlabel("Interaction Type"), Guide.ylabel("Entropy"),
        ColourPalette))

## Here we could plot Entropy and relative rank defficiency

draw(PNG(joinpath("figures", "entropy_v_rank.png"), dpi = 300),
    plot(Outputs,
        x=:RankDefficiencyRel, y=:Entropy,
        color=:InteractionType, alpha = [0.8], Guide.xlabel("Relative rank defficiency"),
        Guide.ylabel("Entropy"),
        ColourPalette))

## Size vs Rank & Entropy

#= TODO
    * Change labels on Y axes (this could come down to changing things in the DataFrame?)
    * Labels on subplots e.g. A/B (see below for an attempt) using layers =#

draw(PNG(joinpath("figures", "size_v_rank&entropy.png"), dpi=300),
    plot(stack(Outputs, [:RankDefficiencyRel, :Entropy], variable_name =:measure, value_name=:value),
        x=:Richness, y=:value,
        color=:InteractionType, ygroup =:measure, Geom.subplot_grid(Geom.point, free_y_axis=true),
        alpha = [0.8], Guide.xlabel("Richness"), Guide.ylabel(nothing),
        ColourPalette))

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

draw(PNG(joinpath("figures", "others_v_entropy.png"), dpi = 300),
    plot(stack(stack(Outputs, [:Entropy], [:Nestedness, :SpectralRadius, :InteractionType],
        variable_name =:measure, value_name=:value), [:Nestedness, :SpectralRadius],
        variable_name =:method, value_name=:metric),
        x=:metric, y=:value,
        color=:InteractionType, ygroup =:measure, xgroup =:method, Geom.subplot_grid(Geom.point, free_y_axis=true),
        alpha = [0.8], Guide.xlabel(nothing), Guide.ylabel(nothing),
        ColourPalette))

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
draw(PNG(joinpath("figures", "entropy_v_AUCall.png"), dpi = 300),
        plot(AUC,
            x =:value, y =:Entropy,
            color=:InteractionType, ygroup =:Type, xgroup =:Dimension,
            Geom.subplot_grid(Geom.point, free_y_axis=true),
            alpha = [0.8], Guide.xlabel("Resilience"), Guide.ylabel("Extinction mechanism"),
            ColourPalette))

## End of script
