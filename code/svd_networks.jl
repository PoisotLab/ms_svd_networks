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
Left as an exercise to the reader 😉
"""
function svd_entropy(N::T) where {T <: DeterministicNetwork}
    F = svd(N)
    Λ = F.S[1:rank(N)]
    λ = Λ ./ sum(Λ)
    return -sum(λ .* log.(λ)) * 1 / log(length(λ))
end

function maxrank(N::T) where {T <: BipartiteNetwork}
    return minimum(size(N))
end

function maxrank(N::T) where {T <: UnipartiteNetwork}
    return richness(N)
end

##FIGURE 1##
###📊 Size vs Rank

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

##💀INCORPORATING EXTINCTIONS

# FUNCTION FROM POISOT 2019
#= This simulates extinctions by removing a RANDOM indiv =#

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

#= modification of extinction function to remove ROW species with the highest number of interactions
    currently this is based on ProbabilityWeights
    ❗TO CHECK: Does it seem 'fair' to use # of interactions or should we look at something like degree()
                Do we want species removal to be proababilistic or 'absolute'
                Maybe need to only make probabilistic when more than one max/min value =#

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

#= we can 'flip' this to remove the LEAST connected species
    here we simply inverse the ProbabilityWeight =#
function extinctions_systematic_least(N::T) where {T <: AbstractBipartiteNetwork}
    # We start by making a copy of the network to extinguish
    Y = [copy(N)];
    # While there is at least one species remaining...
    while richness(last(Y)) > 1
        # sum the number of interactions by col (i.e ⬆ connectance), inverse and convert to ProbabilityWeight
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

#= 📊 SVD Entropy vs prop of spp removed
    (❗need to change the x axis to be more intuitive)
    only using 10 networks for ease of visuals =#
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

###Resilliance - Calculating the area under the curve###
# look at trapezoidrule...
