export SuperCellSystem, Potential, PairPotential, ThreeBodyPotential,
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
#TODO USE ATOMSBASE
#COMMENT OUT UNIT CELL SYSTEM?

struct SuperCellSystem{D,L}
    atoms::StructArray{Atom}
    box_sizes_SC::AbstractVector{L}
end

function SuperCellSystem(positions::AbstractVector{AbstractVector}, masses::AbstractVector,
    box_sizes::AbstractVector, charges = zeros(length(masses))*u"q")
    D = length(box_sizes)
    atoms = StructArray{Atom}((position = positions, mass = masses, charge = charges))
    return SuperCellSystem{D, eltype(box_sizes)}(atoms, box_sizes)
end

function SuperCellSystem(crystal::Crystal{D}) where D
    positions = SimpleCrystals.position(crystal)
    masses = SimpleCrystals.atomic_mass(crystal)
    charges = getindex.(crystal.atoms, :charge)
    atoms = StructArray{Atom}((position = positions,
                               mass = masses,
                               charge = charges))
    box_sizes = norm.(bounding_box(crystal))
    return SuperCellSystem{D, eltype(box_sizes)}(atoms, box_sizes)
end


function SuperCellSystem(df::DataFrame, masses::AbstractVector, box_sizes::AbstractVector,
         x_col, y_col, z_col; charges = zeros(length(masses))*u"q")
    D = length(box_sizes)
    tmp = hcat(df[!,x_col], df[!,y_col], df[!, z_col])
    r = [tmp[i,:] for i in range(1, size(tmp)[1])]
    atoms = StructArray{Atom}((position = r, mass = masses, charge = charges))
    return SuperCellSystem{D, eltype(box_sizes)}(atoms, box_sizes)
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

# struct ForceConstants{V,N,U,T}
#     values::Array{V,N}
#     units::U
#     tol::T
# end

struct FC_val{V,N}
    val::V
    idxs::SVector{Int32,N}
end

function FC_val(val, idxs::Vararg{Integer})
    return FC_val{typeof(val),length(idxs)}(val,SVector(idxs))
end

# struct SparseForceConstants{V,N,U,T}
#     values::AbstractVector{FC_val{V,N}}
#     units::U
#     tol::T
# end

# #type aliases
# const SecondOrderForceConstants{V,U,T} = ForceConstants{V,2,U,T}
# const ThirdOrderForceConstants{V,U,T} = ForceConstants{V,3,U,T}
# const FourthOrderForceConstants{V,U,T} = ForceConstants{V,4,U,T}

# #Construct SparseForceConstants object from dense ForceConstants object 
# function SparseForceConstants(fc::ForceConstants{V,N,U,T}) where {V,N,U,T}
#     #TODO
#     return SparseForceConstants{V,N,U,T}(values, fc.units, fc.tol)
# end