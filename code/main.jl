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
#using Plots
using StatsBase
using Gadfly
using Cairo
using Fontconfig
using DataFrames

## Import the functions and methods we need
include(joinpath(pwd(), "code", "lib", "main.jl"))
include(joinpath(pwd(), "code", "lib", "extinctions.jl"))

## Using the Web Of Life dataset
raw_wol = [web_of_life(x.ID) for x in web_of_life()]
Bs = convert.(BipartiteNetwork, raw_wol)
filter!(B -> richness(B) < 200, Bs)
Ns = convert.(UnipartiteNetwork, Bs)

## 📊 Size vs Rank

#=Plots.plot(scatter(richness.(Bs), rank.(Bs),
        xlabel="Richness", ylabel="Rank", legend=false, dpi=300),
    scatter(richness.(Bs), svd_entropy.(Bs),
        xlabel="Richness", ylabel="Entropy", legend=false, dpi=300))

savefig(joinpath("figures", "size_v_rank.png"))=#

Outputs = DataFrame(Richness = richness.(Bs),
                    Rank = rank.(Bs),
                    Entropy = svd_entropy.(Bs),
                    RankDefficiency = (maxrank.(Bs) .- rank.(Bs)),
                    RankDefficiencyRel = ((maxrank.(Bs) .- rank.(Bs)) ./ maxrank.(Bs)),
                    Nestedness = η.(Bs),
                    SpectralRadiance = ρ.(Bs),
                    InteractionType = [y[:Type_of_interactions] for y in web_of_life()]);

draw(PNG(joinpath("figures", "size_v_rank.png"), 30cm, 10cm, dpi=300),
        hstack(plot(Outputs, x=:Richness, y=:Rank, color=:InteractionType, Geom.point,
            Theme(key_position=:none), alpha = [0.6]),
            plot(Outputs, x=:Entropy, y=:Rank, color=:InteractionType, Geom.point,
            Theme(key_position=:none), alpha = [0.6]),
            plot(Outputs, color=:InteractionType, Geom.blank,
                Guide.colorkey(title="Interaction Type", pos = [0.2,1.5]))))

# 📊 Here we could plot Entropy and various measures of Rank

draw(PNG(joinpath("figures", "entropy_v_rank.png"), 20cm, 20cm, dpi = 300),
    vstack(hstack(plot(Outputs, x=:Richness, y=:Entropy, Geom.point),
                plot(Outputs, x=:Rank, y=:Entropy, Geom.point)),
            hstack(plot(Outputs, x=:RankDefficiency, y=:Entropy, Geom.point, Guide.xlabel("Rank Defficiency")),
                plot(Outputs, x=:RankDefficiencyRel, y=:Entropy, Geom.point, Guide.xlabel("Rank Defficiency (relative)")))))

#=Plots.plot(
    Plots.scatter(richness.(Bs) , svd_entropy.(Bs),
        xlabel="Richness", ylabel="Entropy", legend=false, dpi=300),
    Plots.scatter(rank.(Bs), svd_entropy.(Bs),
        xlabel="Rank", legend=false, dpi=300),
    Plots.scatter(maxrank.(Bs) .- rank.(Bs) , svd_entropy.(Bs),
        xlabel="Rank defficiency", ylabel="Entropy", legend=false, dpi=300),
    Plots.scatter((maxrank.(Bs) .- rank.(Bs)) ./ maxrank.(Bs) , svd_entropy.(Bs),
        xlabel="Rank defficiency (relative)", legend=false, dpi=300)
)

savefig(joinpath("figures", "entropy_v_rank.png"))=#

##📊 COMPARING ENTROPY TO OTHER MEASURES##

draw(PNG(joinpath("figures", "entropy_v_others.png"), 20cm, 20cm, dpi = 300),
    vstack(hstack(plot(Outputs, x=:Nestedness, y=:Entropy, color=:InteractionType, Geom.point,
                    Theme(key_position=:none), alpha = [0.6]),
                plot(Outputs, x=:Nestedness, y=:RankDefficiencyRel, color=:InteractionType, Geom.point,
                    Theme(key_position=:none), alpha = [0.6], Guide.xlabel("Number of Links"))),
            hstack(plot(Outputs, x=:SpectralRadiance, y=:Entropy, color=:InteractionType, Geom.point,
                    Theme(key_position=:none), alpha = [0.6]),
                plot(Outputs, x=:SpectralRadiance, y=:RankDefficiencyRel, color=:InteractionType, Geom.point,
                    Theme(key_position=:none), alpha = [0.6]))))


#=plot(
    # 🐣 Nestedness
    scatter(η.(Bs), svd_entropy.(Bs),
        xlabel="Nestedness", ylabel="Entropy", legend=false, dpi=300),
    # Number of links
    scatter(links.(Bs),svd_entropy.(Bs),
        xlabel="Number of links", legend=false),
    # Connectance
    scatter(connectance.(Bs),svd_entropy.(Bs),
        xlabel="Connectance", legend=false, dpi=300),
    # Linkage density
    scatter(linkage_density.(Bs),svd_entropy.(Bs),
        xlabel="Linkage density", ylabel="Entropy", legend=false, dpi=300),
    # Spectral radius
    scatter(ρ.(Bs), svd_entropy.(Bs),
        xlabel="Spectral radius", ylabel="Entropy", legend=false, dpi=300),
    # Ive added the last two just because I could - doesn't mean I should
    scatter(heterogeneity.(Ns), svd_entropy.(Bs),
        xlabel="Heterogeneity", legend=false, dpi=300),
    scatter(symmetry.(Ns), svd_entropy.(Bs),
        xlabel="Symmetry", legend=false, dpi=300)
)

savefig(joinpath("figures", "entropy_v_others.png"))=#

## 📊 SVD Entropy vs proportion of spp removed

#=plot1 = scatter(title="Most connected")
plot2 = scatter(title="Least connected")
plot3 = scatter(title="Random")
for B in Bs
    Plots.plot!(plot1, (richness(B) .- richness.(extinctions_systematic(B)))/richness(B),
        svd_entropy.(extinctions_systematic(B)), legend=false, c=:grey, alpha=0.5, dpi = 300)
    Plots.plot!(plot2, (richness(B) .- richness.(extinctions_systematic_least(B)))/richness(B),
        svd_entropy.(extinctions_systematic_least(B)), legend=false, c=:grey, alpha=0.5, dpi = 300)
    Plots.plot!(plot3, (richness(B) .- richness.(extinctions(B)))/richness(B), svd_entropy.(extinctions(B)), legend=false, c=:grey, alpha=0.5, dpi = 300)
end

# Remember that plots are objects that can be modified!

for p in [plot1, plot2, plot3]
    yaxis!(p, (0.6, 1.0), "Entropy")
    xaxis!(p, "Proportion of species removed")
end

Plots.plot(plot1, plot2, plot3)

savefig(joinpath("figures", "extinctions_raw.png"))=#

## 📊 Extinction Curves

#=plot4 = scatter(title="Most connected")
plot5 = scatter(title="Least connected")
plot6 = scatter(title="Random")
for B in Bs
    Plots.plot!(plot4, richness.(extinctions_systematic(B))/richness(B),
    (richness(B) .- richness.(extinctions_systematic(B)))/richness(B), legend=false, c=:grey, alpha=0.5, dpi = 300)
    Plots.plot!(plot5, richness.(extinctions_systematic_least(B))/richness(B),
    (richness(B) .- richness.(extinctions_systematic_least(B)))/richness(B), legend=false, c=:grey, alpha=0.5, dpi = 300)
    Plots.plot!(plot6, richness.(extinctions(B))/richness(B),
    (richness(B) .- richness.(extinctions(B)))/richness(B), legend=false, c=:grey, alpha=0.5, dpi = 300)
end
for p in [plot4, plot5, plot6]
    yaxis!(p, "Proprtion of species remaining")
    xaxis!(p, "Proprtion of species removed")
end

for B in Bs
    Plots.plot!(plot4, richness.(extinctions_systematic(B))/richness(B),
    (richness(B) .- richness.(extinctions_systematic(B)))/richness(B), legend=false, c=:grey, alpha=0.5, dpi = 300)
    Plots.plot!(plot5, richness.(extinctions_systematic_least(B))/richness(B),
    (richness(B) .- richness.(extinctions_systematic_least(B)))/richness(B), legend=false, c=:grey, alpha=0.5, dpi = 300)
    Plots.plot!(plot6, richness.(extinctions(B))/richness(B),
    (richness(B) .- richness.(extinctions(B)))/richness(B), legend=false, c=:grey, alpha=0.5, dpi = 300)
end

Plots.plot(plot4, plot5, plot6)

savefig(joinpath("figures", "extinction_curves.png"))=#


## 📊 Resilience vs Rank & Entropy

#=AUC = DataFrame(AUCrandom = [auc((richness(B) .- richness.(extinctions(B)))/richness(B), richness.(extinctions(B))/richness(B)) for B in Bs],
                AUCleast = [auc((richness(B) .- richness.(extinctions_systematic_least(B)))/richness(B), richness.(extinctions_systematic_least(B))/richness(B)) for B in Bs],
                AUCmost = [auc((richness(B) .- richness.(extinctions_systematic(B)))/richness(B), richness.(extinctions_systematic(B))/richness(B)) for B in Bs],
                Rank = ((maxrank.(Bs) .- rank.(Bs)) ./ maxrank.(Bs)),
                Entropy = svd_entropy.(Bs),
                InteractionType = Outputs.InteractionType)
AUC = stack(AUC, [:Rank, :Entropy], variable_name =:measure, value_name=:value)
stack(AUC, [:AUCrandom, :AUCleast, :AUCmost], variable_name =:method, value_name=:auc)

hstack(plot(x = [auc((richness(B) .- richness.(extinctions(B)))/richness(B), richness.(extinctions_systematic(B))/richness(B)) for B in Bs],
            y = svd_entropy.(Bs), color=Outputs.InteractionType, Geom.point,
            Theme(key_position=:none), alpha = [0.6], Guide.xlabel(nothing), Guide.ylabel(nothing), Guide.title("Random")),
        plot(x = [auc((richness(B) .- richness.(extinctions(B)))/richness(B), richness.(extinctions_systematic(B))/richness(B)) for B in Bs],
            y = ((maxrank.(Bs) .- rank.(Bs)) ./ maxrank.(Bs)), color=Outputs.InteractionType, Geom.point,
            Theme(key_position=:none), alpha = [0.6], Guide.xlabel(nothing), Guide.ylabel(nothing)))

plot(AUC, x =:AUCrandom, y =:value,
    color=:InteractionType, xgroup="variable",  Geom.subplot_grid(Geom.point, free_y_axis=true),
    Theme(key_position=:none), alpha = [0.6], Guide.xlabel(nothing), Guide.ylabel(nothing))


plot7 = scatter(title="Most connected")
plot8 = scatter(title="Least connected")
plot9 = scatter(title="Random")
    scatter!(plot7, auc_most.(Bs),svd_entropy.(Bs), legend=false, c=:grey, alpha=0.5, dpi = 300)
    scatter!(plot8, auc_least.(Bs),svd_entropy.(Bs), legend=false, c=:grey, alpha=0.5, dpi = 300)
    scatter!(plot9, auc_random.(Bs),svd_entropy.(Bs), legend=false, c=:grey, alpha=0.5, dpi = 300)
for p in [plot7, plot8, plot9]
    yaxis!(p, "Entropy")
    xaxis!(p, "AUC ≈ Resilience")
end
Plots.plot(plot7, plot8, plot9)

savefig(joinpath("figures", "rank&entropy_v_AUC.png"))=#

## End of script
