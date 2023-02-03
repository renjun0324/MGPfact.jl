module CellTrek
    
using Mamba
using LinearAlgebra: I, Hermitian
using KernelFunctions: kernelmatrix, Kernel
using Distributions: Normal, Beta, truncated
# using KernelFunctions
# using Distributed
# using DataFrames, CSV #just for debuging
# using MLKernels, SpecialFunctions
# using Distributed
# using Distances, KernelFunctions, Mamba
# using Distributions, LinearAlgebra
# using Distributions: rand, length, insupport, logpdf, pdf

include("pseudot.jl")
include("track.jl")
include("modMGPpseudoT.jl")
include("celltracking.jl")

end 