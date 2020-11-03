function optim!(N; f=svd_entropy, target=0.0)
    λ = 0.999
    t0 = 2.0
    Δbest = (f(N)-target)^2.0
    for step in 1:10000
        t = λ^(step-1)*t0
        from = rand(filter(i -> N.A[i], eachindex(N.A)))
        to = rand(filter(i -> !N.A[i], eachindex(N.A)))
        N.A[from], N.A[to] = N.A[to], N.A[from]
        Δstate = (f(N)-target)^2.0
        Δ = Δstate - Δbest
        P = exp(-Δ/t)
        if (rand() ≤ P) & (all(values(degree(N)) .> 0.0))
            Δbest = Δstate
        else
            N.A[from], N.A[to] = N.A[to], N.A[from]
        end
    end
end