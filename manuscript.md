---
bibliography: [references.bib]
---

# Introduction

* The idea of entropy as 'internal complexity' and rank as 'external complexity'
* Why SVD entropy vs 'normal' entropy??
* Relationship between complexity and resilience


This is a citation: @Ginebreda2019QuaEco


# Methods

1. Networks used
 + From web of Life
 + Only bipartite networks (brief description?)
 + removed those with a richness > 200
2. SVD, Entropy in particular
 + normalised

$$ \overline{\lambda_{i}}=\frac{\lambda_{i}}{\Sigma_{i}\lambda_{i}}$$
 + SVD entropy (Shannon type entropy)

$$Entropy = -\frac{1}{\ln(n)}\Big\sum_{i=1}^n \overline{\lambda_{i}}\cdot\ln(\overline{\lambda_{i}})$$
3. Rank
 + Using relative rank deficiency

$$Relative rank deficiency = \frac{Rank_{max} - Rank}{Rank_{max}}$$
4. Other measures of networks
 + Compared to both entropy and relative rank deficiency
 + Something about why these measures
    + Nestedness ($\eta$)
    + Spectral radiance ($\rho$)
5. Extinctions and resilience
 + Species were removed from either the entire network or only along a single dimension i.e. a specific subset of the interacting species such as hosts or pollinators
 + Extinctions were calculated based on three 'mechanisms' by removing; 1) a random individual, 2) the most connected individual and 3) removing the least connected individual.
 + Networks were then simplified after each subsequent species removal by removing all species that no longer had any interactions until all species were removed
 + Comparing the proportion of species remaining to the proportion removed we constructed an extinction curve. Using the area under this curve following the trapezoidal rule we calculates the resilience of the network [ref?]
6. Make note of Julia packages used
 + Is there a ~~lazy~~ smart way to do this from the manifest??
 + link to actual code?

# Results


![This is the legend of the figure](figures/entropy_v_rank.png){#fig:entropy_v_rank}

We can refer to +@fig:entropy_v_rank.

# Discussion

Yes

~~~ julia
for i in eachindex(x)
  x[i] = zero(eltype(x)) # Don't do that
end
~~~

# References
