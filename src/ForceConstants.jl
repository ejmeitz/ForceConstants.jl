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
include("helper_funcs.jl")
include("interactions/LennardJones.jl")
include("interactions/SW.jl")
include("second_order.jl")
include("third_order.jl")
include("fourth_order.jl")
include("dispersion.jl")
include("io.jl")
include("mcc.jl")

include("finite_diff_ifc.jl")
include("autodiff_test.jl")

end
