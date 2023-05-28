export LJ

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

function LJ(σ, ϵ, r_cut; shift_potential = true)
    #Values for U_cut and F_cut are unknown at this point
    temp = LJ{typeof(σ),typeof(ϵ),typeof(r_cut),Integer,Integer}(σ, ϵ, r_cut, 0, 0)

    U_cut = potential(temp, r_cut)
    F_cut = force(temp, r_cut)
    return LJ{typeof(σ),typeof(ϵ),typeof(r_cut),typeof(ϵ),typeof(F_cut)}(σ, ϵ, r_cut, U_cut, F_cut)
end

function potential(pot::LJ, r)
    k = (pot.σ/r)^6
    return 4*pot.ϵ*k*(k-1)
end

function potential_first_deriv(pot::LJ, r)
    return -1*pot.ϵ*(48*(pot.σ^12/r^13) - 24*(pot.σ^6/r^7))
end

function force(pot::LJ, r)
    return -1*potential_first_deriv(pot, r)
end

function potential_second_deriv(pot::LJ, r)
    return 4*pot.ϵ*((12*13*(pot.σ^12/r^14)) - (6*7*(pot.σ^6/r^8)))
end