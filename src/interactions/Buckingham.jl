

struct Buckingham{D,E,F,RC,UC,FC} <: PairPotential
    A::D
    b::L
    c::G
    r_cut::L
    U_cut::UC
    F_cut::FC
end

function Buckingham(A, b, c, r_cut)
    
    @assert unit(1/b) == unit(r_cut) "Units of b and r_cut are not comensurate"

    #Values for U_cut and F_cut are unknown at this point
    temp = Buckingham{typeof(A),typeof(b), typeof(c),typeof(r_cut),Integer,Integer}(
        r_cut, 0, 0)
    
    U_cut = potential(temp, r_cut)
    F_cut = force(temp, r_cut)

    return Buckingham{typeof(A),typeof(b), typeof(c),typeof(r_cut),typeof(U_cut),typeof(F_cut)}(
        r_cut, U_cut, F_cut)
end

energy_unit(pot::Buckingham) = unit(pot.A)
length_unit(pot::Buckingham) = unit(pot.r_cut)


function potential(pot::Buckingham, r)
    return pot.A*exp(-pot.b*r) - pot.c/r^6
end

function potential_nounits(pot::Buckingham, r)
    return ustrip(pot.A*exp(-pot.b*r) - pot.c/r^6)
end

function force(pot::Buckingham, r)
    return pot.A*pot.b*exp(-pot.b*r) - 6*pot.c/r^7
end

function potential_first_deriv(pot::Buckingham, r)
    return 6*pot.c/r^7 - pot.A*pot.b*exp(-pot.b*r)
end

function potential_second_deriv(pot::Buckingham, r)
    return pot.A*pot.b^2*exp(-pot.b*r) - 42*pot.c/r^8
end

function potential_third_deriv(pot::Buckingham, r)
    return 336*pot.c/r^9 - pot.A*pot.b^3*exp(-pot.b*r) 
end

function potential_fourth_deriv(pot::Buckingham, r)
    return pot.A*pot.b^4*exp(-pot.b*r) - 3024*pot.c/r^10 
end