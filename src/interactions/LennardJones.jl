export LJ, energy_unit, length_unit, potential, force

#################
# Lennard Jones #
#################
struct LJ{S,E,R,U,F} <: PairPotential
    σ::S
    ϵ::E
    r_cut::R
    U_cut::U
    F_cut::F
end

function LJ(σ, ϵ, r_cut)
    #Values for U_cut and F_cut are unknown at this point
    temp = LJ{typeof(σ),typeof(ϵ),typeof(r_cut),Integer,Integer}(
        σ, ϵ, r_cut, 0, 0)

    U_cut = potential(temp, r_cut)
    F_cut = force(temp, r_cut)
    return LJ{typeof(σ),typeof(ϵ),typeof(r_cut),typeof(U_cut),typeof(F_cut)}(
        σ, ϵ, r_cut, U_cut, F_cut)
end

energy_unit(pot::LJ) = unit(pot.ϵ)
length_unit(pot::LJ) = unit(pot.σ)

function potential(pot::LJ, r)
    k = (pot.σ/r)^6
    return 4*pot.ϵ*k*(k-1)
end

function potential_nounits(pot::LJ, r)
    k = (ustrip(pot.σ)/r)^6
    return 4*ustrip(pot.ϵ)*k*(k-1)
end

function potential_first_deriv(pot::LJ, r)
    return -1*pot.ϵ*(48*(pot.σ^12/r^13) - 24*(pot.σ^6/r^7))
end

function force(pot::LJ, r)
    return pot.ϵ*(48*(pot.σ^12/r^13) - 24*(pot.σ^6/r^7))
end

function potential_second_deriv(pot::LJ, r)
    return 4*pot.ϵ*((12*13*(pot.σ^12/r^14)) - (6*7*(pot.σ^6/r^8)))
end

function potential_third_deriv(pot::LJ, r)
    return -4*pot.ϵ*((12*13*14(pot.σ^12/r^15)) - (6*7*8(pot.σ^6/r^9)))
end

function potential_fourth_deriv(pot::LJ, r)
    return 4*pot.ϵ*((12*13*14*15(pot.σ^12/r^16)) - (6*7*8*9(pot.σ^6/r^10)))
end