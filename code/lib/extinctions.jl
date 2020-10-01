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
