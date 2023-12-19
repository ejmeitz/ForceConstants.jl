module ForceConstants


using SimpleCrystals
using StructArrays
using Unitful
using StaticArrays
using DataFrames
using DelimitedFiles
using LinearAlgebra
using Combinatorics
using JLD2
using TensorOperations
using Combinatorics
using LoopVectorization
using FastDifferentiation
using cuTENSOR
using CUDA
import CUDA: i32 #TODO: move to weakdep https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)


include("types.jl")
include("interactions/LennardJones.jl")
include("interactions/SW.jl")
include("helper_funcs.jl")
include("third_order.jl")
include("fourth_order.jl")
include("dispersion.jl")
include("io.jl")

# include("finite_diff_ifc.jl")

include("./autodiff/second_order_AD.jl")
include("./autodiff/third_order_AD.jl")
include("./autodiff/fourth_order_AD.jl")

include("./finitediff/second_order_FD.jl")
include("./finitediff/third_order_FD.jl")
include("./finitediff/fourth_order_FD.jl")

include("./analytical/second_order_analytical.jl")
include("./analytical/third_order_analytical.jl")
include("./analytical/fourth_order_analytical.jl")

include("dynmat.jl")
include("mcc.jl")


end
