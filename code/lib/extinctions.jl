"""
    TODO
"""

function extinction(N::T, f::F; dims::Int64=1) where {T <: AbstractBipartiteNetwork, F <: Function}
    if !(dims ∈ [1,2])
        throw(ArgumentError("dims must be 1 or 2 (you used $(dims))"))
    end
    network_series = Vector{T}(undef, richness(N; dims=dims)+1)
    network_series[1] = copy(N)
    extinction_sequence = f(N; dims=dims)
    for (i,sp_to_remove) in enumerate(extinction_sequence)
        species_to_keep = filter(sp -> sp != sp_to_remove, species(network_series[i]; dims=dims))
        K = copy(network_series[i])
        if dims == 1
            K = K[species_to_keep, :]
        else
            K = K[:, species_to_keep]
        end
        network_series[i+1] = simplify(K)
    end
    return network_series
end

"""
    TODO
"""
function _order_species_for_removal(increasing::Nothing=nothing)
    function f(N::T; dims::Int64=1) where {T <: AbstractBipartiteNetwork}
        return StatsBase.shuffle(species(N; dims=dims))
    end
    return f
end

"""
    TODO
"""
function _order_species_for_removal(increasing::Bool=true)
    function f(N::T; dims::Int64=1) where {T <: AbstractBipartiteNetwork}
        k = degree(N; dims=dims)
        kvp = collect(k)
        StatsBase.shuffle!(kvp)
        seq = sort(kvp; by = x -> x[2], rev=!increasing)
        return [s.first for s in seq]
    end
    return f
end

#=

f_random = _order_species_for_removal()
f_increasing = _order_species_for_removal(true)
f_decreasing = _order_species_for_removal(false)

extinction(B, f_random; dims=1)
=#
"""
    TODO
"""
function extinction_robustness(Ns::Vector{T}; dims::Union{Nothing,Int64}=nothing) where {T <: AbstractBipartiteNetwork}
    x = collect(LinRange(0.0, 1.0, length(Ns)))
    if isnothing(dims)
        y = richness.(Ns)./richness(first(Ns))
    else
        y = richness.(Ns; dims=dims)./richness(first(Ns); dims=dims)
    end
    return auc(x, y)
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
