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

    freqs_sq, _ = get_modes(dynmat_sc)
    idxs = sortperm(abs.(freqs_sq))
    freqs_sq[idxs[1:3]] .= 0.0
    ω_THz_sc = sqrt.(freqs_sq .* 4184 .* 1000 .* 1e20)./(2*pi*1e12)

    ### Test Unitcell System ###
    fcc_crystal_1UC = FCC(5.2468u"Å", :Ar, SVector(1,1,1))
    sys_uc = UnitCellSystem(fcc_crystal_1UC, SVector(4,4,4))
    dynmat_uc = dynamicalMatrix(sys_uc, pot, SVector(0.0, 0.0, 0.0)u"Å^-1", 1e-12);

    freqs_sq, _ = get_modes(dynmat_sc)    
    idxs = sortperm(abs.(freqs_sq))
    freqs_sq[idxs[1:3]] .= 0.0
    ω_THz_uc = sqrt.(freqs_sq .* 4184 .* 1000 .* 1e20)./(2*pi*1e12)

    ω_THz_sc = round.(ω_THz_sc, sigdigits = 6)
    ω_THz_uc = round.(ω_THz_uc, sigdigits = 6)
    
    @test all(ω_THz_uc .∈ ω_THz_sc)
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

