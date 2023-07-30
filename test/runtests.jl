using Test
using ForceConstants
using SimpleCrystals
using StaticArrays
using LinearAlgebra


# write tests here

@testset "Supercell vs Unitcell" begin
    
    pot = LJ(3.4u"Å", 0.24037u"kcal * mol^-1", 8.5u"Å")

    ### Test Supercell System ###
    fcc_crystal = FCC(5.2468u"Å", :Ar, SVector(4,4,4))
    sys_sc = SuperCellSystem(fcc_crystal)
    dynmat_sc = dynamicalMatrix(sys_sc, pot, 1e-12);

    freqs_sq, _ = get_modes(dynmat_sc, 3)

    ω_THz_sc = sqrt.(freqs_sq .* 4184 .* 1000 .* 1e20)./(2*pi*1e12)

    ### Test Unitcell System ###
    fcc_crystal_1UC = FCC(5.2468u"Å", :Ar, SVector(1,1,1))
    sys_uc = UnitCellSystem(fcc_crystal_1UC, SVector(4,4,4))
    dynmat_uc = dynamicalMatrix(sys_uc, pot, SVector(0.0, 0.0, 0.0)u"Å^-1", 1e-12);

    freqs_sq, _ = get_modes(dynmat_uc, 3)    
    ω_THz_uc = sqrt.(freqs_sq .* 4184 .* 1000 .* 1e20)./(2*pi*1e12)

    ω_THz_sc = round.(ω_THz_sc, sigdigits = 6)
    ω_THz_uc = round.(ω_THz_uc, sigdigits = 6)
    
    @test issubset(ω_THz_uc, ω_THz_sc)
end

@testset "MCC3 Blocked" begin
    
    #Code assumes first element is the one that got passed in
    @test collect(multiset_permutations([1,2,3],3))[1] == [1,2,3]
    @test collect(multiset_permutations([1,2,2],3))[1] == [1,2,2]
    @test collect(multiset_permutations([:i,:k,:i],3))[1] == [:i,:k,:i]

    #Load test dataset
    phi, dynmat, F3, K3_actual, freqs_sq = load("./test_data/TEP.jld2", "phi", "dynmat", "F3", "K3", "freqs_sq")
    N_modes = length(freqs_sq) #should be 96

    n_blocks = 3
    K3_full = mcc3(CuArray(F3), CuArray(phi))
    K3_blocked = mcc3(CuArray(F3), CuArray(phi),n_blocks)

    @test isapprox(K3_actual, K3_full, atol = 1e-6)
    @test isapprox(K3_actual, K3_blocked, atol = 1e-6)
    
end


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

