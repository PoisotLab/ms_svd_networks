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
include(joinpath(pwd(), "code", "lib", "annealing.jl"))

S = 14
minL = S+4
maxL = S^2-4
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

relcon = (N) -> (links(N).-EcologicalNetworks.richness(N; dims=1))./(prod(size(N)).-EcologicalNetworks.richness(N; dims=1))
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

p1 = Gadfly.plot(sims, x=:relcon, y=:entropy, color=:optim, Stat.smooth(smoothing=0.5), Geom.line, Scale.color_discrete_manual(colorant"black", colorant"grey"), Guide.colorkey(title="Type", labels=["Most complex","Least complex"]))
push!(p1, layer(sims, x=:relcon, y=:entropy, color=:optim, alpha=[0.2]))
push!(p1, Coord.cartesian(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0))
push!(p1, Guide.xlabel("Connectance"))
push!(p1, Guide.ylabel("SVD entropy"))

p2 = Gadfly.plot(sims, x=:relcon, y=:reldef, color=:optim, Stat.smooth(smoothing=0.2), Geom.line, Scale.color_discrete_manual(colorant"black", colorant"grey"), Guide.colorkey(title="Type", labels=["Most complex","Least complex"]))
push!(p2, layer(sims, x=:relcon, y=:reldef, color=:optim, alpha=[0.1]))
push!(p2, Coord.cartesian(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0))
push!(p2, Guide.xlabel("Connectance"))
push!(p2, Guide.ylabel("Rel. rank def."))

Gadfly.set_default_plot_size(21cm, 12cm)
draw(
    PNG(joinpath("figures", "minmax_combined.png"), dpi=300),
    hstack(p1, p2)
)
Gadfly.set_default_plot_size(18cm, 12cm)
