import LinearAlgebra
import EcologicalNetworks
import Statistics
import Plots
import StatsBase # import instead of use

#Using the Web Of Life dataset
ids = [x.ID for x in web_of_life()]
Bs = convert.(BipartiteNetwork, [web_of_life(id) for id in ids])
filter!(B -> richness(B) < 200, Bs)

bipart = Bs
unipart = convert.(UnipartiteNetwork, bipart)

function LinearAlgebra.rank(N::T) where {T<:DeterministicNetwork}
    return rank(N.A)
end

function LinearAlgebra.svd(N::T) where {T<:DeterministicNetwork}
    return svd(N.A)
end

# Entropy based on SVD
function svd_entropy(N::T) where {T<:DeterministicNetwork}
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

plot(scatter(richness.(bipart), rank.(bipart),
        xlabel = "Richness", ylabel = "Rank", legend = false, dpi = 300),
    scatter(richness.(bipart), svd_entropy.(bipart),
        xlabel = "Richness", ylabel = "Entropy", legend = false, dpi = 300))

#📊 Here we could plot Entropy and various measures of Rank
plot(
    scatter(richness.(bipart) , svd_entropy.(bipart),
        xlabel = "Richness", ylabel = "Entropy", legend = false, dpi = 300),
    scatter(rank.(bipart), svd_entropy.(bipart),
        xlabel = "Rank", legend = false, dpi = 300),
    scatter(maxrank.(bipart) .- rank.(bipart) , svd_entropy.(bipart),
        xlabel = "Rank defficiency", ylabel = "Entropy", legend = false, dpi = 300),
    scatter((maxrank.(bipart) .- rank.(bipart))./maxrank.(bipart) , svd_entropy.(bipart),
        xlabel = "Rank defficiency (relative)", legend = false, dpi = 300)
)

##📊 COMPARING ENTROPY TO OTHER MEASURES##

plot(
    #🐣 Nestedness
    scatter(η.(bipart), svd_entropy.(bipart),
        xlabel = "Nestedness", ylabel = "Entropy", legend = false, dpi = 300),
    #Number of links
    scatter(links.(bipart),svd_entropy.(bipart),
        xlabel = "Number of links", legend = false),
    #Connectance
    scatter(connectance.(bipart),svd_entropy.(bipart),
        xlabel = "Connectance", legend = false, dpi = 300),
    #Linkage density
    scatter(linkage_density.(bipart),svd_entropy.(bipart),
        xlabel = "Linkage density", ylabel = "Entropy", legend = false, dpi = 300),
    #Ive added the last two just because I could - doesn't mean I should
    scatter(heterogeneity.(unipart), svd_entropy.(bipart),
        xlabel = "Heterogeneity", legend = false, dpi = 300),
    scatter(symmetry.(unipart), svd_entropy.(bipart),
        xlabel = "Symmetry", legend = false, dpi = 300)
)

##💀INCORPORATING EXTINCTIONS

#FUNCTION FROM POISOT 2019
#= This simulates extinctions by removing a RANDOM indiv =#

function extinctions(N::T) where {T<:AbstractBipartiteNetwork}
    # We start by making a copy of the network to extinguish
    Y = [copy(N)]
    # While there is at least one species remaining...
    while richness(last(Y)) > 1
        # We remove one species randomly
        remain = sample(
            species(last(Y); dims = 2),
            richness(last(Y); dims = 2) - 1,
            replace = false,
        )
        # Remaining species
        R = last(Y)[:, remain]
        simplify!(R)
        # Then add the simplified network (without the extinct species) to our collection
        push!(Y, copy(R))
    end
    return Y
end

#=modification of extinction function to remove ROW species with the highest number of interactions
    currently this is based on ProbabilityWeights
    ❗TO CHECK: Does it seem 'fair' to use # of interactions or should we look at something like degree()
                Do we want species removal to be proababilistic or 'absolute'
                Maybe need to only make probabilistic when more than one max/min value=#

function extinctions_systematic(N::T) where {T<:AbstractBipartiteNetwork}
    # We start by making a copy of the network to extinguish
    Y = [copy(N)];
    # While there is at least one species remaining...
    while richness(last(Y)) > 1
        #sum the number of interactions by row (i.e ⬆ connectance) and convert to ProbabilityWeight
        w = ProbabilityWeights(vec(sum(last(Y), dims = 1)./sum(sum(last(Y), dims = 1))))
        # We remove the species based on probability weighting (w)
        remain =  sample(
            species(last(Y); dims = 2),
            w,
            richness(last(Y); dims = 2) - 1,
            replace = false,
        )
        # Remaining species
        R = last(Y)[:, remain]
        simplify!(R)
        # Then add the simplified network (without the extinct species) to our collection
        push!(Y, copy(R))
    end
    return Y
end

#=we can 'flip' this to remove the LEAST connected species
    here we simply inverse the ProbabilityWeight=#
function extinctions_systematic_least(N::T) where {T<:AbstractBipartiteNetwork}
    # We start by making a copy of the network to extinguish
    Y = [copy(N)];
    # While there is at least one species remaining...
    while richness(last(Y)) > 1
        #sum the number of interactions by col (i.e ⬆ connectance), inverse and convert to ProbabilityWeight
        w = ProbabilityWeights(vec(1 ./ sum(last(Y), dims = 1)./sum(sum(last(Y), dims = 1))))
        # We remove the species based on probability weighting (w)
        remain =  sample(
            species(last(Y); dims = 2),
            w,
            richness(last(Y); dims = 2) - 1,
            replace = false,
        )
        # Remaining species
        R = last(Y)[:, remain]
        simplify!(R)
        push!(Y, copy(R))
    end
    return Y
end

#=📊 SVD Entropy vs prop of spp removed
    (❗need to change the x axis to be more intuitive)
    only using 10 networks for ease of visuals =#
scatter()
for i = 1:10
    scatter!((richness.(extinctions_systematic(bipart[i]); dims = 2) / richness(bipart[i]; dims = 2)),
        svd_entropy.(extinctions_systematic(bipart[i])), legend = false)
end
plot3 = scatter!(title = "Systematic - most", xlabel = "Propotion of species removed")
scatter()
for i = 1:10
    scatter!(richness.(extinctions_systematic_least(bipart[i])) / richness(bipart[i]),
        svd_entropy.(extinctions_systematic_least(bipart[i])), legend = false)
end
plot2 = scatter!(title = "Systematic - least", ylabel = "Entropy")
scatter()
for i = 1:10
    scatter!(richness.(extinctions(bipart[i])) / richness(bipart[i]), svd_entropy.(extinctions(bipart[i])), legend = false)
end
plot1 = scatter!(title = "Random")

plot(plot1, plot2, plot3)

###Resilliance - Calculating the area under the curve###
#look at trapezoidrule...
