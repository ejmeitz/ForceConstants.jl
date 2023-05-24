struct Atom{P,C,M}
    position::P
    mass::M
    charge::C
end

function Atom(position, mass; charge = 0.0u"q")
    return Atom{typeof(position), typeof(mass), typeof(charge)}(position, mass, charge)
end

#####################################################

abstract type System end

struct UnitCellSystem{D,L} <: System
    atoms::StructArray{Atom}
    num_unit_cells::SVector{D,Int}
    L::L
end

struct SuperCellSystem{D,L} <: System
    atoms::StructArray{Atom}
    L::L
end

masses(sys::System) = sys.atoms.mass
mass(sys::System, i::Int) = sys.atoms.mass[i]
positions(sys::System) = sys.atoms.position
position(sys::System, i ::Int) = sys.atoms.position[i]
charges(sys::System) = sys.atoms.charge
charge(sys::System) = sys.atoms.charge[i]
n_atoms(sys::System) = length(sys.atoms)
#######################################################

abstract type Potential end
abstract type PairPotential <: Potential end

########################################################

