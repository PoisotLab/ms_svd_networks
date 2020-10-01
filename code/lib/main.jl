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
