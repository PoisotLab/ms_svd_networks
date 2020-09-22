import LinearAlgebra
import EcologicalNetworks
import Statistics
import Plots
import StatsBase # import instead of use, this will make the global namespace cleaner

#Using the Web Of Life dataset
ids = [x.ID for x in web_of_life()]
Bs = convert.(BipartiteNetwork, [web_of_life(id) for id in ids])
filter!(B -> richness(B) < 200, Bs)

bipart = Bs

#= A very powerful design pattern in Julia is "overloading",
where we create new methods for our own types -
so to make the code easier, we will define a method for svd and rank!=#

function LinearAlgebra.rank(N::T) where {T<:DeterministicNetwork}
    return rank(N.A)
end

function LinearAlgebra.svd(N::T) where {T<:DeterministicNetwork}
    return svd(N.A)
end


##FIGURE 1##
#Size vs Rank

scatter(richness.(bipart), rank.(bipart), leg = false, aspectratio = 1);
xaxis!("Richness", (0, 200));
yaxis!("Rank", (0, 100));


# Entropy based on SVD
function svd_entropy(N::T) where {T<:DeterministicNetwork}
    F = svd(N)
    Λ = F.S[1:rank(N)]
    λ = Λ ./ sum(Λ)
    return -sum(λ .* log.(λ)) * 1 / log(length(λ))
end

#= NOTE If you want to calculate more things, you should do more functions,
each of which does one thing -
or you can look into the CSV and DataFrame paackages,
and use this as a way to store your results.=#

#Here we could plot Rank and Entropy
scatter(rank.(bipart), svd_entropy.(bipart))
xlabel!("Rank");
ylabel!("Entropy");

#= NOTE This is also a variant of having a specific method depending on the type of arguments
 - a.k.a. dispatch - it's a super important mechanism in Julia
 (and multiple dispatch even more so!) =#

function maxrank(N::T) where {T <: BipartiteNetwork}
    return minimum(size(N))
end

function maxrank(N::T) where {T <: UnipartiteNetwork}
    return richness(N)
end

scatter(maxrank.(bipart) .- rank.(bipart) , svd_entropy.(bipart))
xlabel!("Rank defficiency")
ylabel!("Entropy")


scatter((maxrank.(bipart) .- rank.(bipart))./maxrank.(bipart) , svd_entropy.(bipart))
xlabel!("Rank defficiency (relative)")
ylabel!("Entropy")

##COMPARING ENTROPY TO OTHER MEASURES##

plot(
    #🐣 Nestedness
    scatter(η.(bipart), svd_entropy.(bipart),
        xlabel = "Nestedness", ylabel = "Entropy", legend = false),
    #↕ Number of links
    scatter(links.(bipart),svd_entropy.(bipart),
        xlabel = "Number of links", legend = false),
    #🔀 Connectance
    scatter(connectance.(bipart),svd_entropy.(bipart),
        xlabel = "Connectance", ylabel = "Entropy", legend = false),
    #⇵ Linkage density
    scatter(linkage_density.(bipart),svd_entropy.(bipart),
        xlabel = "Linkage density", legend = false)
)

##INCORPORATING EXTINCTIONS

#FUNCTION FROM POISOT 2019
#= This simulates extinctions
by removing a RANDOM indiv =#

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

scatter()
for i = 1:length(bipart)
    scatter!(richness.(extinctions(bipart[i])) / richness(bipart[i]), svd_entropy.(extinctions(bipart[i])), legend = false)
end
scatter!()
xlabel!("Proportion of species removed")
ylabel!("Entropy")

#=modification to remove ROW species with the highest number of interactions
    currently this is based on ProbabilityWeights
    ❗TO CHECK: Does it seem 'fair' to use # of interactions or should we look at something like degree()
                Do we want species removal to be proababilistic or 'absolute'=#
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
        #sum the number of interactions by row (i.e ⬆ connectance) and convert to ProbabilityWeight
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
        # Then add the simplified network (without the extinct species) to our collection
        push!(Y, copy(R))
    end
    return Y
end

scatter()
for i = 1:length(bipart)
    scatter!(richness.(extinctions_systematic(bipart[i])) / richness(bipart[i]), svd_entropy.(extinctions_systematic(bipart[i])), legend = false)
end
scatter!()
xlabel!("Proportion of species removed")
ylabel!("Entropy")

scatter()
for i = 1:length(bipart)
    scatter!(richness.(extinctions_systematic_least(bipart[i])) / richness(bipart[i]), svd_entropy.(extinctions_systematic_least(bipart[i])), legend = false)
end
scatter!()
xlabel!("Proportion of species removed")
ylabel!("Entropy")
