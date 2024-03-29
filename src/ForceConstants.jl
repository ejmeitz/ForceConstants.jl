module ForceConstants

using SimpleCrystals
using StructArrays
using Unitful
using StaticArrays
using DataFrames #& only import whats needed?
using DelimitedFiles
using LinearAlgebra
using Combinatorics
using JLD2
using TensorOperations
using Combinatorics

#*only needed by AD parts
using FastDifferentiation
using Symbolics
using RuntimeGeneratedFunctions

using cuTENSOR
using CUDA
import CUDA: i32 #TODO: move to weakdep https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)


include("types.jl")
include("interactions/LennardJones.jl")
include("interactions/SW.jl")

include("helper_funcs.jl")
include("asr.jl")

sw_generated_dir = joinpath(@__DIR__, "autodiff", "generated_derivatives", "SW")
include.(filter(contains(r".jl$"), readdir(sw_generated_dir; join=true)))
# include("./autodiff/autodiff_helper.jl")
include("./autodiff/second_order_AD.jl")
include("./autodiff/third_order_AD.jl")
# include("./autodiff/fourth_order_AD.jl")

include("./finitediff/finitediff_helper.jl")
include("./finitediff/second_order_FD.jl")
include("./finitediff/third_order_FD.jl")
# include("./finitediff/fourth_order_FD.jl")
# include("./finitediff/check_ifc.jl")

include("./analytical/second_order_analytical.jl")
include("./analytical/third_order_analytical.jl")
# include("./analytical/fourth_order_analytical.jl")

# include("gruneisen.jl")
include("ifc.jl")
include("dynmat.jl")
include("mcc.jl")
include("io.jl") 


end
