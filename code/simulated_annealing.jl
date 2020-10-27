import Pkg
Pkg.activate("code")

using LinearAlgebra, Statistics
using EcologicalNetworks
using StatsBase
using Gadfly
import Cairo, Fontconfig
using Colors
using DataFrames

include(joinpath(pwd(), "code", "lib", "main.jl"))
include(joinpath(pwd(), "code", "lib", "extinctions.jl"))

f = (x) -> x ./ sum(x)

# Examples
function optim!(N; f=svd_entropy, target=0.0)
    λ = 0.999
    t0 = 2.0
    #target = 0.0
    Δbest = (svd_entropy(N)-target)^2.0
    #out = [(Δbest, t0)]
    for step in 1:10000
        t = λ^(step-1)*t0
        from = rand(filter(i -> N.A[i], eachindex(N.A)))
        to = rand(filter(i -> !N.A[i], eachindex(N.A)))
        N.A[from], N.A[to] = N.A[to], N.A[from]
        Δstate = (svd_entropy(N)-target)^2.0
        Δ = Δstate - Δbest
        P = exp(-Δ/t)
        if (rand() ≤ P) & (all(values(degree(N)) .> 0.0))
            Δbest = Δstate
        else
            N.A[from], N.A[to] = N.A[to], N.A[from]
        end
        #push!(out, (Δbest, t))
    end
end


S = 12

minL = S+2
maxL = S^2
nreps = 20

L = minL:2:maxL
Ns = []
for l in L
    m = zeros(Bool, (S, S))
    m[sample(eachindex(m), l, replace=false)] .= true
    N = BipartiteNetwork(m)
    while !(all(values(degree(N)) .> 0.0))
        m = zeros(Bool, (S, S))
        m[sample(eachindex(m), l, replace=false)] .= true
        N = BipartiteNetwork(m)
    end
    push!(Ns, N)
end

Ns = repeat(Ns, nreps)

MAX = copy.(Ns)
MIN = copy.(Ns)

optim!.(MAX; target=1.0)
optim!.(MIN; target=0.0)

sims = DataFrame(
    relcon = Float64[],
    optim = Symbol[],
    entropy = Float64[],
    reldef = Float64[]
)

relcon = (N) -> (links(N).-richness(N; dims=1))./(prod(size(N)).-richness(N; dims=1))
reldef = (N) -> ((maxrank(N) - rank(N)) / maxrank(N))

for i in 1:length(MAX)
    push!(sims, (
        relcon(MAX[i]),
        :maximum,
        svd_entropy(MAX[i]),
        reldef(MAX[i])
        )
    )
    push!(sims, (
        relcon(MIN[i]),
        :minimum,
        svd_entropy(MIN[i]),
        reldef(MIN[i])
        )
    )
end

p1 = plot(sims, x=:relcon, y=:entropy, color=:optim, Stat.smooth, Geom.line, Scale.color_discrete_hue())
push!(p1, layer(sims, x=:relcon, y=:entropy, color=:optim, alpha=[0.2]))
push!(p1, Coord.cartesian(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0))

p2 = plot(sims, x=:relcon, y=:reldef, color=:optim, Stat.smooth, Geom.line, Scale.color_discrete_hue())
push!(p2, layer(sims, x=:relcon, y=:reldef, color=:optim, alpha=[0.2]))
push!(p2, Coord.cartesian(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0))

hstack(p1, p2)