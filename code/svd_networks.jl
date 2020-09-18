using LinearAlgebra
using EcologicalNetworks
using Statistics
using Plots
import Mangal # import instead of use, this will make the global namespace cleaner

# NOTE I have changed the code to use the web of life dataset
# This removes the need to use Mangal for development
ids = [x.ID for x in web_of_life()]
Bs = convert.(BipartiteNetwork, [web_of_life(id) for id in ids])
filter!(B -> richness(B) < 200, Bs)

bipart = Bs

# A very powerful design pattern in Julia is "overloading", where we create new methods for our own types - so to make the code easier, we will define a method for svd and rank!

function LinearAlgebra.rank(N::T) where {T<:DeterministicNetwork}
    return rank(N.A)
end

function LinearAlgebra.svd(N::T) where {T<:DeterministicNetwork}
    return svd(N.A)
end


##FIGURE 1##
#Size vs Rank

scatter(richness.(bipart), rank.(bipart), leg = false, aspectratio = 1)
xaxis!("Richness", (0, 200))
yaxis!("Rank", (0, 100))


# Entropy based on SVD
# NOTE this function does one thing, and one thing only: it returns a value for the entropy of singular values up to the rank of the matrix
function svd_entropy(N::T) where {T<:DeterministicNetwork}
    F = svd(N)
    Λ = F.S[1:rank(N)]
    λ = Λ ./ sum(Λ)
    return -sum(λ .* log.(λ)) * 1 / log(length(λ))
end

#CURRENTLY THIS OUTPUT IS A MESS...
# NOTE If you want to calculate more things, you should do more functions, each of which does one thing - or you can look into the CSV and DataFrame paackages, and use this as a way to store your results.

#Here we could plot Rank and Entropy
scatter(rank.(bipart), svd_entropy.(bipart))
xlabel!("Rank");
ylabel!("Entropy");

# NOTE This is also a variant of having a specific method depending on the type of arguments - a.k.a. dispatch - it's a super important mechanism in Julia (and multiple dispatch even more so!)

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

#=
1 Extract and convert to bipartate network
2 Calcluate nestedness - save ouptu
3 Connectance? - save output
=#


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

#Caluclaute richness of an extinction matrix
function extinction_richness(N::T) where {T<:BinaryNetwork}
    return richness(N; dims = 1)
end



plot()
for i = 1:5
    extinct = extinctions(bipart[i])
    Extinction = (Prop = Any[], Entropy = Any[])
    for j = 1:length(extinct)
        #decompose matrix
        F = svd(extinct[j].A)
        λ = F.S[1:rank(extinct[j].A)] / sum(F.S[1:rank(extinct[j].A)])
        push!(Extinction.Entropy, -sum(λ_entropy.(λ)))
        push!(Extinction.Prop, (extinction_richness(extinct[j]) / richness(bipart)))
    end
    plot!(Extinction.Prop, Extinction.Entropy, legend = false)
end
plot!()
