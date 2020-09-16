using LinearAlgebra
using EcologicalNetworks
using Statistics
using Plots
using Mangal

#Function to 'normalise'
#=this would calulate the 'contibution to total entopy'
of each σ value =#
function λ_entropy(λ)
    λ*log2(λ)
end

#pull 100 datasets from Mangal
eco_networks = networks();

unipart = convert.(UnipartiteNetwork, eco_networks[1:5]#=only using five for testing=#);
#= 🚨 Keeps throwing error when running full collection of eco_networks
'has both in and out degree'
❗TO DO -> add a checkpoint to see if criterea are met =#
bipart = convert.(BipartiteNetwork, unipart);


##FIGURE 1##
#Size vs Rank

Rank = Any[];
Rich = Any[];

for i in 1:length(bipart)
    #❗TO CHECK: should we be working with the unipart/bipartite network
    push!(Rich, richness(bipart[i]));
    push!(Rank, rank(unipart[i].A));
end

plot(Rich, Rank, seriestype = :scatter, legend = false);
xlabel!("Rank");
ylabel!("Species Richness");


##A SIMPLIFIED APPROACH TO SVD AND ENTROPY CALCULATIONS##

function SVD_Entropy(N::T) where {T <: UnipartiteNetwork}

    bipart = convert(BipartiteNetwork, N);
    #decompose matrix
    F = svd(bipart.A);
    #normalise σ values -> ❗TO CHECK: HOW WE USE RANK (of raw matrix/SVD)
    λ = F.S[1:rank(bipart.A)]/sum(F.S[1:rank(bipart.A)]);
    #Calculate total Entropy by apply function along the vector
    Entropy =  -sum(λ_entropy.(λ));
    #Calculate Rank differential
    RankDiff = length(F.S) - rank(bipart.A);
    #Record Original Rank
    Rank = rank(bipart.A);
    #Return outputs ->
    #❗TO DO: how to make the outputs 'nicer'/intuitive
    return (λ, Entropy, RankDiff, Rank)
end


#CURRENTLY THIS OUTPUT IS A MESS...
#Calc entropy etc for the various networks
X = SVD_Entropy.(unipart);

#=
#Here we could plot Rank and Entropy
plot(X[1][2],
    X[1][2],
    title = "", dpi=1300);
xlabel!("Rank diff");
ylabel!("Entropy");=#

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

function extinctions(N::T) where {T <: AbstractBipartiteNetwork}
# We start by making a copy of the network to extinguish
Y = [copy(N)]
# While there is at least one species remaining...
while richness(last(Y)) > 1
    # We remove one species randomly
    remain = sample(species(last(Y); dims=2), richness(last(Y); dims=2)-1, replace=false)
    # Remaining species
    R = last(Y)[:,remain]
    simplify!(R)
    # Then add the simplified network (without the extinct species) to our collection
    push!(Y, copy(R))
end
return Y
end

#Caluclaute richness of an extinction matrix
function extinction_richness(N::T) where
{T<:BinaryNetwork}
return richness(N; dims=1)
end



plot()
for i in 1:5
    extinct  = extinctions(bipart[i]);
    Extinction = (Prop = Any[],Entropy = Any[]);
    for j in 1:length(extinct)
        #decompose matrix
        F = svd(extinct[j].A);
        λ = F.S[1:rank(extinct[j].A)]/sum(F.S[1:rank(extinct[j].A)]);
        push!(Extinction.Entropy, -sum(λ_entropy.(λ)));
        push!(Extinction.Prop, (extinction_richness(extinct[j])/richness(bipart)));
    end
    plot!(Extinction.Prop,Extinction.Entropy, legend = false);
end
plot!()
