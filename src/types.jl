export UnitCellSystem, SuperCellSystem,
    masses, mass, position, positions, charges, charge, n_atoms, n_atoms_per_uc

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
    atoms::Dict{Tuple,StructArray{Atom}}
    num_unit_cells::SVector{D,Int}
    box_sizes_UC::L
    box_sizes_SC::L
end

function UnitCellSystem(crystal::Crystal{D}, num_unit_cells::SVector{D,Int}) where D
    @assert sum(crystal.N_unit_cells) == D "Crystal must be a single unit cell"

    unit_cell_indices = Iterators.product([1:num_unit_cells[i] for i in range(1,D)]...)

    positions = SimpleCrystals.position(crystal)
    masses = SimpleCrystals.atomic_mass(crystal)
    charges = getindex.(crystal.atoms, :charge)
    box_sizes = norm.(eachrow(bounding_box(crystal)))
    box_sizes_SC = norm.(num_unit_cells .* eachrow(bounding_box(crystal)))
    #Generate data for all atoms in all unit cells
    atoms = Dict()

    for uc_idx in unit_cell_indices
        positions_new = []
        unitcell_origin = [box_sizes[1]*(uc_idx[1]-1), box_sizes[2]*(uc_idx[2]-1), box_sizes[3]*(uc_idx[3]-1)]
        for i in range(1, length(crystal))
            push!(positions_new, positions[i] .+ unitcell_origin)
        end
        atoms[uc_idx] = StructArray{Atom}((position = positions_new, mass = masses, charge = charges))
    end

    return UnitCellSystem{D, typeof(box_sizes)}(atoms, num_unit_cells, box_sizes, box_sizes_SC)
end

mass(sys::UnitCellSystem, atom_idx::Int) = sys.atoms[(1,1,1)].mass[atom_idx]
position(sys::UnitCellSystem, uc_idx::Tuple, atom_idx::Int) = sys.atoms[uc_idx].position[atom_idx]
charge(sys::UnitCellSystem, atom_idx::Int) = sys.atoms[(1,1,1)].charge[atom_idx]
n_atoms(sys::UnitCellSystem) = length(sys.atoms)*length(sys.atoms[(1,1,1)])
n_atoms_per_uc(sys::UnitCellSystem) = length(sys.atoms[(1,1,1)])

struct SuperCellSystem{D,L} <: System
    atoms::StructArray{Atom}
    box_sizes_SC::L
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

function SuperCellSystem(r0, masses, box_sizes::AbstractVector; charges = zeros(length(masses))*u"q")
    D = length(box_sizes)
    atoms = StructArray{Atom}((position = r0, mass = masses, charge = charges))
    return SuperCellSystem{D, typeof(box_sizes)}(atoms, box_sizes)
end

function SuperCellSystem(df::DataFrame, masses::AbstractVector, box_sizes::AbstractVector,
         x_col, y_col, z_col; charges = zeros(length(masses))*u"q")
    D = length(box_sizes)
    tmp = hcat(df[!,x_col], df[!,y_col], df[!, z_col])
    r = [tmp[i,:] for i in range(1, size(tmp)[1])]
    atoms = StructArray{Atom}((position = r, mass = masses, charge = charges))
    return SuperCellSystem{D, typeof(box_sizes)}(atoms, box_sizes)
end

atomkeys(sys::SuperCellSystem) = keys(sys.atoms[1])
hasatomkey(sys::SuperCellSystem, x::Symbol) = x ∈ atomkeys(sys)
masses(sys::SuperCellSystem) = sys.atoms.mass
mass(sys::SuperCellSystem, i::Int) = sys.atoms.mass[i]
positions(sys::SuperCellSystem) = sys.atoms.position
position(sys::SuperCellSystem, i::Int) = sys.atoms.position[i]
positions_1D(sys::SuperCellSystem) = vcat(positions(sys))
charges(sys::SuperCellSystem) = hasatomkey(sys, :charge) ? sys.atoms.charge : throw(ArgumentError("Charge is not a key"))
charge(sys::SuperCellSystem, i::Int) =  hasatomkey(sys, :charge) ? sys.atoms.charge[i] : throw(ArgumentError("Charge is not a key"))
n_atoms(sys::SuperCellSystem) = length(sys.atoms)
#######################################################

abstract type Potential end
abstract type PairPotential <: Potential end
abstract type ThreeBodyPotential <: Potential end

########################################################

