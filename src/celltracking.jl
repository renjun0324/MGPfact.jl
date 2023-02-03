function celltracking(model, yx, Q, L, rt, iterations, nc)
    P = size(yx)[1]
    SC = Dict{Symbol, Any}(
        # :Yx => Matrix(yx["cl"]'),
        # :N => N,
        :Yx => Matrix(yx'),
        :P => P,
        :Q => Q,
        :L => L,
        :rt => rt,
        :sigma => 1,
        :delta => 1
    )   
    println(1)

    inits = [
        Dict{Symbol, Any}(
            :Yx => SC[:Yx][1:Q,:],
            :T => rand(SC[:P]),
            :X => rand(truncated(Normal(0, 1), -1.0, 1.0), SC[:P]),
            :rho => 0.5,
            :p => [rand(Beta(1, 1)) for i in 1:SC[:L], t in 1:SC[:P]],
            :Tb => zeros(SC[:L]) .+ 0.2,
            :s2_t => 1,
            :s2_x => 1,
            :m_t => 0.5,
            :lambda1 => rand(truncated(Normal(0, 99), 0.0, 1.0), SC[:L]),
            :lambda2 => rand(truncated(Normal(0, 99), 0.0, 0.000001), SC[:L]),
            :lambda3 => rand(truncated(Normal(0, 99), 0.0, 1.0), SC[:L])        
        )
    for i in 1:nc
    ]

    scheme = [AMWG(:T, 10),
              Slice(:Tb, 1, transform=true),
              Slice(:X, 1, transform=true),
              Slice(:p, 1, transform=true),
              Slice(:rho, 1, transform=true),
              Slice(:m_t, 1),
              Slice([:s2_t, :s2_x], [100, 100]),
              Slice(:lambda1, 1, transform=true),
              Slice(:lambda2, 1, transform=true),
              Slice(:lambda3, 1, transform=true)]

    # setinputs!(model, SC)
    # setinits!(model, inits)
    # setsamplers!(model, scheme)
    # model = CellTrek.modMGPpseudoT()
    @time sim1 = mcmc(model, SC, inits, iterations, burnin = 0, chains = nc)
    return sim1
end