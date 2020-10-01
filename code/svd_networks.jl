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


# For this type of projets, calling packages with using is fine, as the
# namespace is never going to get very big - one situation where you would want
# to use import instead is if you were using a package only in a single line
# (for example, using DataFrames but import CSV is fairly common)

## Using the Web Of Life dataset
raw_wol = [web_of_life(x.ID) for x in web_of_life()]
Bs = convert.(BipartiteNetwork, raw_wol)
filter!(B -> richness(B) < 200, Bs)
Ns = convert.(UnipartiteNetwork, Bs)

## Redefine methods for rank and svd

"""
    LinearAlgebra.rank(N::T) where {T <: DeterministicNetwork}

Returns the rank of the adjacency matrix of a deterministic ecological network.
"""
function LinearAlgebra.rank(N::T) where {T <: DeterministicNetwork}
    return rank(N.A)
end

"""
    LinearAlgebra.svd(N::T) where {T <: DeterministicNetwork}

Performs Singular Value Decomposition on the adjacency matrix of a deterministic network
"""
function LinearAlgebra.svd(N::T) where {T <: DeterministicNetwork}
    return svd(N.A)
end

## Entropy based on SVD

"""
    svd_entropy(N::T) where {T <: DeterministicNetwork}

Returns the entropy of the adjacency matrix of a deterministic network using the singular values from a Singualr Value Decomposition
"""
function svd_entropy(N::T) where {T <: DeterministicNetwork}
    F = svd(N)
    Λ = F.S[1:rank(N)]
    λ = Λ ./ sum(Λ)
    return -sum(λ .* log.(λ)) * 1 / log(length(λ))
end

"""
    maxrank(N::T) where {T <: BipartiteNetwork}

Returns the maximum possible rank of a Bipartite Network
"""
function maxrank(N::T) where {T <: BipartiteNetwork}
    return minimum(size(N))
end

"""
    maxrank(N::T) where {T <: UnipartiteNetwork}

Returns the maximum possible rank of a Unipartite Network
"""
function maxrank(N::T) where {T <: UnipartiteNetwork}
    return richness(N)
end

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

## 💀 INCORPORATING EXTINCTIONS

# FUNCTION FROM POISOT 2019
"""
    extinctions(N::T) where {T <: AbstractBipartiteNetwork}

This returns a vector of bipartite ecological networks simulating an extinction trajectory, where each network is the result of
removing a random individual from the preceeding network
"""

function extinctions(N::T) where {T <: AbstractBipartiteNetwork}
    # We start by making a copy of the network to extinguish
    Y = [copy(N)]
    # While there is at least one species remaining...
    while richness(last(Y)) > 1
        # We remove one species randomly
        remain = sample(
            species(last(Y); dims=2),
            richness(last(Y); dims=2) - 1,
            replace=false,
        )
        # Remaining species
        R = last(Y)[:, remain]
        simplify!(R)
        # Then add the simplified network (without the extinct species) to our collection
        push!(Y, copy(R))
    end
    return Y
end

"""
    extinctions_systematic(N::T) where {T <: AbstractBipartiteNetwork}

This returns a vector of bipartite ecological networks simulating an extinction trajectory, where each network is the result of
removing the most connected individual from the preceeding network
"""
function extinctions_systematic(N::T) where {T <: AbstractBipartiteNetwork}
    # We start by making a copy of the network to extinguish
    Y = [copy(N)];
    # While there is at least one species remaining...
    while richness(last(Y)) > 1
        # sum the number of interactions by row (i.e ⬆ connectance) and convert to ProbabilityWeight
        w = ProbabilityWeights(vec(sum(last(Y), dims=1) ./ sum(sum(last(Y), dims=1))))
        # We remove the species based on probability weighting (w)
        remain =  sample(
            species(last(Y); dims=2),
            w,
            richness(last(Y); dims=2) - 1,
            replace=false,
        )
        # Remaining species
        R = last(Y)[:, remain]
        simplify!(R)
        # Then add the simplified network (without the extinct species) to our collection
        push!(Y, copy(R))
    end
    return Y
end

"""
    extinctions_systematic_least(N::T) where {T <: AbstractBipartiteNetwork}

This returns a vector of bipartite ecological networks simulating an extinction trajectory, where each network is the result of
removing the least connected individual from the preceeding network
"""
function extinctions_systematic_least(N::T) where {T <: AbstractBipartiteNetwork}
    # We start by making a copy of the network to extinguish
    Y = [copy(N)];
    # While there is at least one species remaining...
    while richness(last(Y)) > 1
        # sum the number of interactions by col (i.e ⬆ connectance), INVERSE and convert to ProbabilityWeight
        w = ProbabilityWeights(vec(1 ./ sum(last(Y), dims=1) ./ sum(sum(last(Y), dims=1))))
        # We remove the species based on probability weighting (w)
        remain =  sample(
            species(last(Y); dims=2),
            w,
            richness(last(Y); dims=2) - 1,
            replace=false,
        )
        # Remaining species
        R = last(Y)[:, remain]
        simplify!(R)
        push!(Y, copy(R))
    end
    return Y
end

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

## Calculating the Area under the species extinction curve (AUC)

#=Following the trapezoidal rule we can construct ⟂ trapezoids between known points of the curve
    and calculate the sum of their areas
    The area of a trapezoid can be defined as:
                    Area = 0.5(Δ𝑥)(𝑦ᵢ+𝑦ᵢ₋₁)
                         where Δ𝑥 = 𝑥ᵢ - 𝑥ᵢ₋₁
    In this instance 𝑦 would be the proportion of species remaining in the network and
    𝑥 the proportion of species remaining in the network =#

#=
    Using this logic we would use the 'inner' values of (for example) the vector of 𝑦 values
    twice and the 'outer' values only once. Thus we could have two copies of the 𝑦-axis vector
    but remove the first element from one vector and the final element from the second vector and
    calculate the sum of the pairs. This would also apply for the 𝑥-axis values. We could then calculate
    the area of the individual trapezoids and then summate said areas to get the total AUC.
=#

#= NOTE from TP

It's a lot better to define one function for the AUC, and then apply it to the
simulation. Writing a general, re-usable function is a huge time saver, and also
a best practice. You should re-write the functions to use this.

=#

"""
    TODO
"""
function auc(x::Vector{T}, y::Vector{T}) where {T <: Number}
    @assert length(x) == length(y)
    area = 0.0
    for i in 2:length(x)
        area += 0.5 * (x[i] - x[i-1])*(y[i]+y[i-1])
    end
    return area
end


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
