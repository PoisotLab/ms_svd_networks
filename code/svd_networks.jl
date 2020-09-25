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
using Plots
using StatsBase


# For this type of projets, calling packages with using is fine, as the
# namespace is never going to get very big - one situation where you would want
# to use import instead is if you were using a package only in a single line
# (for example, using DataFrames but import CSV is fairly common)

## Using the Web Of Life dataset
raw_wol = [web_of_life(x.ID) for x in web_of_life()]
Bs = convert.(BipartiteNetwork, raw_wol)
filter!(B -> richness(B) < 200, Bs)
Ns = convert.(UnipartiteNetwork, Bs)

# This is a convention from EcologicalNetworks, but I mostly tend to use B for
# bipartite, N for unipartite or any general network, R for random, and P for
# probabilist - I have reformated the code so that it follows this convention,
# this is also less typing

## Redefine methods for rank and svd

# The new methods work on the abstract (union? can't remember how it's defined)
# type DeterministicNetwork --- the syntax with """ is a standard way to
# document a function in Julia, so if you do ?rank, you will now see these
# methods

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

plot(scatter(richness.(Bs), rank.(Bs),
        xlabel="Richness", ylabel="Rank", legend=false, dpi=300),
    scatter(richness.(Bs), svd_entropy.(Bs),
        xlabel="Richness", ylabel="Entropy", legend=false, dpi=300))

savefig(joinpath("figures", "size_v_rank.png"))

# 📊 Here we could plot Entropy and various measures of Rank
plot(
    scatter(richness.(Bs) , svd_entropy.(Bs),
        xlabel="Richness", ylabel="Entropy", legend=false, dpi=300),
    scatter(rank.(Bs), svd_entropy.(Bs),
        xlabel="Rank", legend=false, dpi=300),
    scatter(maxrank.(Bs) .- rank.(Bs) , svd_entropy.(Bs),
        xlabel="Rank defficiency", ylabel="Entropy", legend=false, dpi=300),
    scatter((maxrank.(Bs) .- rank.(Bs)) ./ maxrank.(Bs) , svd_entropy.(Bs),
        xlabel="Rank defficiency (relative)", legend=false, dpi=300)
)

savefig(joinpath("figures", "entropy_v_rank.png"))

##📊 COMPARING ENTROPY TO OTHER MEASURES##

plot(
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

savefig(joinpath("figures", "entropy_v_others.png"))

## 💀INCORPORATING EXTINCTIONS

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
    (❗need to change the x axis to be more intuitive) =#
plot1 = scatter(title="Most connected")
plot2 = scatter(title="Least connected")
plot3 = scatter(title="Random")
for B in Bs
    plot!(plot1, (richness.(extinctions_systematic(B); dims=2) / richness(B; dims=2)),
        svd_entropy.(extinctions_systematic(B)), legend=false, c=:grey, alpha=0.5)
    plot!(plot2, richness.(extinctions_systematic_least(B)) / richness(B),
        svd_entropy.(extinctions_systematic_least(B)), legend=false, c=:grey, alpha=0.5)
    plot!(plot3, richness.(extinctions(B)) / richness(B), svd_entropy.(extinctions(B)), legend=false, c=:grey, alpha=0.5)
end

# Remember that plots are objects that can be modified!

for p in [plot1, plot2, plot3]
    yaxis!(p, (0.6, 1.0), "Entropy")
end

plot(plot1, plot2, plot3)

savefig(joinpath("figures", "extinctions_raw.png"))

## 📊 Extinction Curves

plot4 = scatter(title="Most connected")
plot5 = scatter(title="Least connected")
plot6 = scatter(title="Random")
for B in Bs
    plot!(plot4, richness.(extinctions_systematic(B))/richness(B),
    (richness(B) .- richness.(extinctions_systematic(B)))/richness(B), legend=false, c=:grey, alpha=0.5)
    plot!(plot5, richness.(extinctions_systematic_least(B))/richness(B),
    (richness(B) .- richness.(extinctions_systematic_least(B)))/richness(B), legend=false, c=:grey, alpha=0.5)
    plot!(plot6, richness.(extinctions(B))/richness(B),
    (richness(B) .- richness.(extinctions(B)))/richness(B), legend=false, c=:grey, alpha=0.5)
end
for p in [plot4, plot5, plot6]
    yaxis!(p, "Proprtion of species remaining")
    xaxis!(p, "Proprtion of species removed")
end

plot(plot4, plot5, plot6)

savefig(joinpath("figures", "extinction_curves.png"))

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

"""
    auc_most(N::T) where {T <: DeterministicNetwork}

Returns the area under a species extinction curve where the most connected species are removed, this is calcualted following the trapezoidal rule
"""
function auc_most(N::T) where {T <: DeterministicNetwork}
    #Calculate the proprotion of species REMAINING i.e. the 𝑌 axis
    Y = richness.(extinctions_systematic(N))/richness(N)
    #Calculate the proprtion of species REMOVED i.e. the 𝘟 axis
    X = (richness(N) .- richness.(extinctions_systematic(N)))/richness(N)
    return sum((collect(Iterators.rest(Y, 2)) #= 'remove' the first element =# + collect(Iterators.take(Y, length(Y)-1)#= 'remove' the final element =#)).*
    ((collect(Iterators.rest(X, 2)) - collect(Iterators.take(X, length(X)-1)))./2))
end

"""
    auc_least(N::T) where {T <: DeterministicNetwork}

Returns the area under a species extinction curve where the least connected species are removed, this is calcualted following the trapezoidal rule
"""
function auc_least(N::T) where {T <: DeterministicNetwork}
    #Calculate the proprotion of species REMAINING i.e. the 𝑌 axis
    Y = richness.(extinctions_systematic_least(N))/richness(N)
    #Calculate the proprtion of species REMOVED i.e. the 𝘟 axis
    X = (richness(N) .- richness.(extinctions_systematic_least(N)))/richness(N)
    return sum((collect(Iterators.rest(Y, 2)) + collect(Iterators.take(Y, length(Y)-1))).*
    ((collect(Iterators.rest(X, 2)) - collect(Iterators.take(X, length(X)-1)))./2))
end

"""
    auc_random(N::T) where {T <: DeterministicNetwork}

Returns the area under a species extinction curve where random species are removed, this is calcualted following the trapezoidal rule
"""
function auc_random(N::T) where {T <: DeterministicNetwork}
    #Calculate the proprotion of species REMAINING i.e. the 𝑌 axis
    Y = richness.(extinctions(N))/richness(N)
    #Calculate the proprtion of species REMOVED i.e. the 𝘟 axis
    X = (richness(N) .- richness.(extinctions(N)))/richness(N)
    return sum((collect(Iterators.rest(Y, 2)) + collect(Iterators.take(Y, length(Y)-1))).*
    ((collect(Iterators.rest(X, 2)) - collect(Iterators.take(X, length(X)-1)))./2))
end

## 📊 Resilience vs Entropy

plot7 = scatter(title="Most connected")
plot8 = scatter(title="Least connected")
plot9 = scatter(title="Random")
    scatter!(plot7, auc_most.(Bs),svd_entropy.(Bs), legend=false, c=:grey, alpha=0.5)
    scatter!(plot8, auc_least.(Bs),svd_entropy.(Bs), legend=false, c=:grey, alpha=0.5)
    scatter!(plot9, auc_random.(Bs),svd_entropy.(Bs), legend=false, c=:grey, alpha=0.5)
for p in [plot7, plot8, plot9]
    yaxis!(p, "Entropy")
    xaxis!(p, "AUC ≈ Resilience")
end
plot(plot7, plot8, plot9)

savefig(joinpath("figures", "entropy_v_AUC.png"))

## End of script
