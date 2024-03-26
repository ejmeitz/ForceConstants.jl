using Test
using ForceConstants
using SimpleCrystals
using StaticArrays
using LinearAlgebra
using CUDA
using Combinatorics
using JLD2
using Unitful

# Allow CUDA device to be specified
const DEVICE = get(ENV, "DEVICE", "0")
run_gpu_tests = false
if CUDA.functional()
    device!(parse(Int, DEVICE))
    run_gpu_tests = true
    @info "The GPU tests will be run on device $DEVICE"
else
    @warn "The GPU tests will not be run as a CUDA-enabled device is not available"
end

#& Check that calculating sparse 3rd equals generating it from the dense
#& Various U_TEP tests, per mode vs total vs parallel vs serial

# pot = LJ(3.4u"Å", 0.24037u"kcal * mol^-1", 8.5u"Å")
# fcc_crystal = FCC(5.2468u"Å", :Ar, SVector(4,4,4))
# sys = SuperCellSystem(fcc_crystal)
# dynmat = dynamicalMatrix(sys, pot, 1e-12);
# freqs_sq, phi = get_modes(dynmat, 3)
# ifc3 = third_order_IFC(sys, pot, 1e-12)
# m = 39.95
# ifc3_mw = mass_weight_third_order!(ifc3, m*ones(length(fcc_crystal)))
# K3 = mcc3(CuArray{Float32}(ifc3_mw.values), CuArray{Float32}(phi), 256, 1e-12)

# @testset "Supercell vs Unitcell" begin
    
#     pot = LJ(3.4u"Å", 0.24037u"kcal * mol^-1", 8.5u"Å")

#     ### Test Supercell System ###
#     fcc_crystal = FCC(5.2468u"Å", :Ar, SVector(4,4,4))
#     sys_sc = SuperCellSystem(fcc_crystal)
#     dynmat_sc = dynamicalMatrix(sys_sc, pot, 1e-12);

#     freqs_sq, _ = get_modes(dynmat_sc, 3)

#     ω_THz_sc = sqrt.(freqs_sq .* 4184 .* 1000 .* 1e20)./(2*pi*1e12)

#     ### Test Unitcell System ###
#     fcc_crystal_1UC = FCC(5.2468u"Å", :Ar, SVector(1,1,1))
#     sys_uc = UnitCellSystem(fcc_crystal_1UC, SVector(4,4,4))
#     dynmat_uc = dynamicalMatrix(sys_uc, pot, SVector(0.0, 0.0, 0.0)u"Å^-1", 1e-12);

#     freqs_sq, _ = get_modes(dynmat_uc, 3)    
#     ω_THz_uc = sqrt.(freqs_sq .* 4184 .* 1000 .* 1e20)./(2*pi*1e12)

#     ω_THz_sc = round.(ω_THz_sc, sigdigits = 6)
#     ω_THz_uc = round.(ω_THz_uc, sigdigits = 6)
    
#     @test issubset(ω_THz_uc, ω_THz_sc)
# end

@testset "MCC3, MCC3 Blocked" begin
    
    #Load test dataset
    if run_gpu_tests
        m = 39.95
        phi, dynmat, F3, K3_actual, freqs_sq = load("./test_data/TEP.jld2", "phi", "dynmat", "F3", "K3", "freqs_sq")
        N_modes = length(freqs_sq) #should be 96

        block_size = 32
        F3 ./= sqrt(m^3)

        K3_full = Array(mcc3(CuArray{Float32}(F3), CuArray{Float32}(phi)))
        K3_blocked = Array(mcc3(CuArray{Float32}(F3), CuArray{Float32}(phi), block_size))

        @test isapprox(K3_actual, K3_full, atol = 1e-6)
        @test isapprox(K3_actual, K3_blocked, atol = 1e-6)
    end
    
end

@testset "MCC vs MCC Lowmem" begin
    
    #Load test dataset
    if run_gpu_tests
        m = 39.95
        phi, dynmat, F3, K3_actual, freqs_sq = load("./test_data/TEP.jld2", "phi", "dynmat", "F3", "K3", "freqs_sq")
        N_modes = length(freqs_sq) #should be 96

        block_size = 32
        F3 ./= sqrt(m^3)
        cuF3 = CuArray{Float32}(F3)
        cuPhi = CuArray{Float32}(phi)
        K3_full = Array(mcc3_tens_opt(cuF3, cuPhi))
        K3_blocked = Array(mcc3(cuF3, cuPhi, block_size))

        K3_lowmem1 = Array(mcc3(cuF3, cuPhi))
        K3_lowmem2 = Array(mcc3!(cuF3, cuPhi))


        @test isapprox(K3_actual, K3_full, atol = 1e-6)
        @test isapprox(K3_actual, K3_blocked, atol = 1e-6)
        @test isapprox(K3_actual, K3_lowmem1, atol = 1e-6)
        @test isapprox(K3_actual, K3_lowmem2, atol = 1e-6)
    end

end

@testset "LJ (FD, AD, Analytical)" begin
    
    pot = LJ(3.4u"Å", 0.24037u"kcal * mol^-1", 8.5u"Å")
    fcc_crystal = FCC(5.2468u"Å", :Ar, SVector(4,4,4))
    sys = SuperCellSystem(fcc_crystal);
    
    tol = 1e-12
    calc_AD = AutoDiffCalculator(tol, pot.r_cut)
    calc_FD = FiniteDiffCalculator(tol, pot.r_cut)
    calc_analytical = AnalyticalCalculator(tol, pot.r_cut)

    ifc2_AD = second_order(sys, pot, calc_AD)
    ifc2_FD = second_order(sys, pot, calc_FD)
    ifc2_analytical = second_order(sys, pot, calc_analytical)

    @test isapprox(ifc2_AD.values, ifc2_FD.values, rtol = 0.01)
    @test isapprox(ifc2_analytical.values, ifc2_AD.values, atol = 1e-12)

    ifc3_AD = third_order(sys, pot, calc_AD)
    # ifc3_FD = third_order(sys, pot, calc_FD)
    ifc3_analytical = third_order(sys, pot, calc_analytical)

    # @test isapprox(ifc3_AD.values, ifc3_FD.values, rtol = 0.01)
    @test isapprox(ifc3_analytical.values, ifc3_AD.values, atol = 1e-12)

end

# @test "SW" begin

#     # Analytical vs Sparse vs FD vs AD (Pair)
#     # Sparse vs FD vs AD (Three-Body)

# end



# ω_THz = get_dispersion_points(sys_uc, pot; directions = ([1.0, 0.0, 0.0],), tol = 1e-12);


# #Since I hard coded 1,0,0 only plot x part
# #plotting like this is also stupid
# scatter()
# for key in keys(ω_THz)
#     k_vals = key[1] .* ones(length(ω_THz[key]))
#     scatter!(k_vals, ω_THz[key], markersize = 4, mc = :blue, legend = :none) #this is stupid but it works
# end
# title!("FCC Dispersion [1,0,0] Direction")
# xlabel!("k_x")
# ylabel!("ω [THz]")
# savefig("C:/Users/ejmei/Desktop/test.png")


