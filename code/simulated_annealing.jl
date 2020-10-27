import Pkg
Pkg.activate("code")

using LinearAlgebra, Statistics
using EcologicalNetworks
using StatsBase
using Gadfly
import Cairo, Fontconfig
using Colors

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

minL = 12+2
maxL = S^2-2
nreps = 25

L = repeat(minL:2:maxL, nreps)
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

MAX = copy.(Ns)
MIN = copy.(Ns)

optim!.(MAX; target=1.0)
optim!.(MIN; target=0.0)

reldef = (N) -> ((maxrank(N) - rank(N)) / maxrank(N))

ls = unique(L)
m = zeros(Float64, length(ls))
M = zeros(Float64, length(ls))
mr = zeros(Float64, length(ls))
Mr = zeros(Float64, length(ls))

for (i,l) in enumerate(ls)
    idx = findall(L .== l)
    m[i] = first(findmin(svd_entropy.(MIN[idx])))
    M[i] = first(findmax(svd_entropy.(MAX[idx])))
    mr[i] = mean(reldef.(MIN[idx]))
    Mr[i] = mean(reldef.(MAX[idx]))
end


X = (ls.-(S-1))./(S^2-(S-1))

p2 = plot(x=X, y=M, Geom.path)
push!(p2, layer(x=X, y=m, Geom.path))
push!(p2, layer(x=X, ymin=m, ymax=M, Geom.ribbon))
push!(p2, Guide.xlabel("Corrected connectance"))
push!(p2, Guide.ylabel("Entropy"))

draw(
    PNG(joinpath("figures", "minmax_entropy.png"), dpi = 300),
    p2
)


p3 = plot(x=X, y=Mr, Geom.path, color=[colorant"black"])
push!(p3, layer(x=X, y=mr, Geom.path, color=[colorant"grey"]))
push!(p3, Guide.xlabel("Corrected connectance"))
push!(p3, Guide.ylabel("Average rank defficiency"))

draw(
    PNG(joinpath("figures", "minmax_rank.png"), dpi = 300),
    p3
)

draw(
    PNG(joinpath("figures", "minmax_combined.png"), dpi = 300),
    hstack(p2, p3)
)
