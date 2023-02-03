function modMGPpseudoT()
    modMGPpseudoT = Model(
        Yx = Stochastic(2,
            (T, Tb, C, X, lambda, L, P, Q, sigma, delta) -> 
                MultivariateDistribution[(
                        S = [(
                                tb = Tb[i];
                                cl = round.(C[i,:]);
                                k = pseudot.bifurcateKernel(tb, lambda[i,:], sigma, 0.0, delta);
                                x = [pseudot.Trejactory(T[j], cl[j], X[j]) for j in 1:P] ;
                                kernelmatrix(k, x)
                            )
                            for i in 1:L
                            ];     
                        MvNormal(zeros(P), Matrix(Hermitian(sum(S) .+ 1E-6 .* Matrix(I, P, P))))            )
                for q in 1:Q
                ]
        ),
        rho = Stochastic(
            () -> truncated(Beta(1.0, 1.0), 1E-8, 0.99999999),
        ),
        a = Logical(2,
            (T, Tb, L, P) -> [
                2.0 - 2.0 / (1.0 + exp(-1*(T[t]-Tb[i])) )
            for i in 1:L, t in 1:P
            ]
        ),
        p = Stochastic(2,
            (rho, T, Tb, L, P, a) -> 
                begin
                    UnivariateDistribution[(
                        Beta(a[i,t] * (1-rho), a[i,t] * rho))
                    for i in 1:L, t in 1:P
                    ]
                end
        ),
        C = Logical(2,
            (p, L, P, rho) -> [
                argmax([p[i,t], 1-p[i,t]])
            for i in 1:L, t in 1:P
            ]
        ),
        T = Stochastic(1,
            (m_t, s2_t, P, rt) ->
                UnivariateDistribution[
                    # if t == rt
                    #     truncated(Normal(0, 0.001), 0.0, 1.0)
                    # # elseif t == ed
                    # #     truncated(Normal(m_t, sqrt(s2_t)), 0.9, 1.0)
                    # else
                    #     truncated(Normal(m_t, sqrt(s2_t)), 0.0, 1.0)
                    # end
                    if t in rt
                        truncated(Normal(0, 0.01), 0.0, 1.0)
                    else
                        truncated(Normal(m_t, sqrt(s2_t)), 0.0, 1.0)
                    end
                    
                for t in 1:P
                ]
        ),
        X = Stochastic(1,
            (s2_x, P) ->
                UnivariateDistribution[
                    truncated(Normal(0, sqrt(s2_x)), -1.0, 1.0)
                for t in 1:P
                ]
        ),
        Tb = Stochastic(1,
            (L) -> UnivariateDistribution[
                truncated(Gamma(i, 1/L), 0.0, 1.0)
            for i in 1:L
            ]
        ),
        m_t = Stochastic(
            () -> truncated(Normal(0.5, 99), 0, 1)
        ),
        s2_t = Stochastic(
            () -> truncated(InverseGamma(0.01, 0.01), 1.0, Inf)
        ),
        s2_x = Stochastic(
            () -> truncated(InverseGamma(0.01, 0.01), 1.0, Inf)
        ),
        lambda1 = Stochastic(1,
            () -> truncated(Normal(0, 99), 0.0, 1.0)
        ),
        lambda2 = Stochastic(1,
            () -> truncated(Normal(0, 99), 0.0, 0.000001)
        ),
        lambda3 = Stochastic(1,
            () -> truncated(Normal(0, 99), 0.0, 1.0)
        ),
        lambda = Logical(2,
            (lambda1, lambda2, lambda3) ->
                [lambda1 lambda2 lambda3]
        )
    )
    return modMGPpseudoT
end

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