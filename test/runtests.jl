
# using Mamba, RData, Distributions
# nc = 1
# L = 3
# Q = 3
# iterations = 200
# rt = [1,2,3]

# yx = RData.load("test/murp_matrix_pca.rda")
# yx = yx["murp_matrix_pca"][:,1:Q]
# P = size(yx)[1]

# SC = Dict{Symbol, Any}(
#     # :Yx => Matrix(yx["cl"]'),
#     # :N => N,
#     :Yx => Matrix(yx'),
#     :P => P,
#     :Q => Q,
#     :L => L,
#     :rt => rt,
#     :sigma => 1,
#     :delta => 1
# )   

# inits = [
#     Dict{Symbol, Any}(
#         :Yx => SC[:Yx][1:Q,:],
#         :T => rand(SC[:P]),
#         :X => rand(truncated(Normal(0, 1), -1.0, 1.0), SC[:P]),
#         :rho => 0.5,
#         :p => [rand(Beta(1, 1)) for i in 1:SC[:L], t in 1:SC[:P]],
#         :Tb => zeros(SC[:L]) .+ 0.2,
#         :s2_t => 1,
#         :s2_x => 1,
#         :m_t => 0.5,
#         :lambda1 => rand(truncated(Normal(0, 99), 0.0, 1.0), SC[:L]),
#         :lambda2 => rand(truncated(Normal(0, 99), 0.0, 0.000001), SC[:L]),
#         :lambda3 => rand(truncated(Normal(0, 99), 0.0, 1.0), SC[:L])        
#     )
# for i in 1:nc
# ]

# scheme = [AMWG(:T, 10),
#           Slice(:Tb, 1, transform=true),
#           Slice(:X, 1, transform=true),
#           Slice(:p, 1, transform=true),
#           Slice(:rho, 1, transform=true),
#           Slice(:m_t, 1),
#           Slice([:s2_t, :s2_x], [100, 100]),
#           Slice(:lambda1, 1, transform=true),
#           Slice(:lambda2, 1, transform=true),
#           Slice(:lambda3, 1, transform=true)]

# model = CellTrek.modMGPpseudoT()
# setinputs!(model, SC)
# setinits!(model, inits)
# setsamplers!(model, scheme)

# @time sim1 = mcmc(model, SC, inits, iterations, burnin = 0, chains = nc)
# write(string("2_pseudotime/2.1_julia_result/iter",sim1.model.iter,"_bi",sim1.model.burnin,".jls"), sim1)

# using JLD2
# @save string("2_pseudotime/2.1_julia_result/inits.jld2") inits

#------------------------------------------------------------------------------
#
#                                 1 thread
#
#------------------------------------------------------------------------------

nc = 1
L = 3
Q = 3
rt = 1
iterations = 10

using CellTrek, Mamba, RData
test_data = joinpath(dirname(pathof(CellTrek)),"..","test","test.rda")
yx = RData.load(test_data)
yx = yx["murp_matrix_pca"][:,1:Q]

# running model
model = CellTrek.modMGPpseudoT()
SC,inits,scheme = CellTrek.Initialize(yx, Q, L, rt, iterations, nc)
setinputs!(model, SC)
setinits!(model, inits)
setsamplers!(model, scheme)
@time sim = mcmc(model, SC, inits, iterations, burnin = 0, chains = nc)

# save
using JLD2
write(string("iter",sim.model.iter,"_bi",sim.model.burnin,".jls"), sim)
@save string("inits.jld2") inits

#------------------------------------------------------------------------------
#
#                                 mutiple threads
#
#------------------------------------------------------------------------------

# using RData
# nc = 3
# L = 3
# Q = 3
# rt = [1,2,3]
# yx = RData.load("test/murp_matrix_pca.rda")
# yx = yx["murp_matrix_pca"][:,1:Q]
# iterations = 200

# using Distributed
# addprocs(nc)
# @everywhere using CellTrek, Mamba, RData, Distributions
# # sim = CellTrek.celltracking(model, yx, Q, L, rt, 100, nc)
# model = CellTrek.modMGPpseudoT()
# SC,inits,scheme = CellTrek.Initialize(yx, Q, L, rt, iterations, nc)
# setinputs!(model, SC)
# setinits!(model, inits)
# setsamplers!(model, scheme)
# @time sim1 = mcmc(model, SC, inits, iterations, burnin = 0, chains = nc)

# save
# using JLD2
# write(string("iter",sim.model.iter,"_bi",sim.model.burnin,".jls"), sim)
# @save string("inits.jld2") inits
