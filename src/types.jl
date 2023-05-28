export UnitCellSystem, SuperCellSystem,
    masses, mass, position, positions, charges, charge, n_atoms

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
    box_sizes::L
end

function UnitCellSystem(crystal::Crystal)
    positions = SimpleCrystals.position(crystal)
    masses = SimpleCrystals.atomic_mass(crystal)
    charges = getindex.(crystal.atoms, :charge)
    atoms = StructArray{Atom}((position = positions,
                               mass = masses,
                               charge = charges))
    box_sizes = norm.(eachrow(bounding_box(crystal)))
    return SuperCellSystem{D, typeof(box_sizes)}(atoms, crystal.N_unit_cells, box_sizes)
end

struct SuperCellSystem{D,L} <: System
    atoms::StructArray{Atom}
    box_sizes::L
end

function SuperCellSystem(crystal::Crystal{D}) where D
    positions = SimpleCrystals.position(crystal)
    masses = SimpleCrystals.atomic_mass(crystal)
    charges = getindex.(crystal.atoms, :charge)
    atoms = StructArray{Atom}((position = positions,
                               mass = masses,
                               charge = charges))
    box_sizes = norm.(eachrow(bounding_box(crystal)))
    return SuperCellSystem{D, typeof(box_sizes)}(atoms, box_sizes)
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

