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
    tmp = hcat(df[!,x_col], df[!,y_col], df[!, z_col])
    r = [tmp[i,:] for i in range(1, size(tmp)[1])]
    atoms = StructArray{Atom}((position = r, mass = masses, charge = charges))
    return SuperCellSystem{3, eltype(box_sizes)}(atoms, box_sizes)
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

abstract type AbstractForceConstants{O} end

struct DenseForceConstants{O,V,U,T} <: AbstractForceConstants{O}
    values::Array{V,O}
    units::U
    tol::T
end

Base.size(fc::DenseForceConstants) = size(fc.values)
Base.getindex(fc::DenseForceConstants{O}, idxs::Vararg{Integer, O}) where O = fc.values[idxs...]
n_modes(fc::DenseForceConstants) = size(fc)[1]

#type aliases
const SecondOrderForceConstants{V,U,T} = DenseForceConstants{2,V,U,T}
const ThirdOrderForceConstants{V,U,T} = DenseForceConstants{3,V,U,T}
const FourthOrderForceConstants{V,U,T} = DenseForceConstants{4,V,U,T}

struct FC_val{V,O}
    val::V
    idxs::Vector{Int32}
end
value(fcv::FC_val) = fcv.val
idx(fcv::FC_val) = fcv.idxs

function FC_val(val, idxs::Vararg{Integer, O}) where O
    return FC_val{typeof(val), O}(val, Vector{Int32}(idxs))
end


struct SparseForceConstants{O,V,U,T} <: AbstractForceConstants{O}
    values::AbstractVector{FC_val{V,O}} #& MAKE THIS A STRUCT ARRAY???
    units::U
    tol::T
end

# struct SparseForceConstantsTest{O,V,U,T} <: AbstractForceConstant{O}
#     values::StructArray{FC_val{V,O}} #& MAKE THIS A STRUCT ARRAY???
#     units::U
#     tol::T
# end

Base.length(fc::SparseForceConstants) = length(fc.values)

#type aliases
const SparseSecondOrder{V,U,T} = SparseForceConstants{2,V,U,T}
const SparseThirdOrder{V,U,T} = SparseForceConstants{3,V,U,T}
const SparseFourthOrder{V,U,T} = SparseForceConstants{4,V,U,T}


# #Construct SparseForceConstants object from dense ForceConstants object 
function SparseForceConstants(fc::ForceConstants{O,V,U,T}) where {O,V,U,T}
    N_modes = size(fc)[1]
    
    idxs = with_replacement_combinations(1:N_modes, O)
    values = FC_Val{V,O}[]
    for idx in idxs
        if fc[idx...] != 0.0
            for perm in permutations(idx)
                push!(values, FC_Val(fc[perm...], perm...))
            end
        end
    end

    return SparseForceConstants{O,V,U,T}(values, fc.units, fc.tol)
end

# Storage is multiple vectors, acts as single vector
# Allows for multi-threaded construction
struct MultiVectorStorage{T}
    data::Vector{Vector{T}}
end

function MultiVectorStorage(dt::DataType, N_vecs)
    return MultiVectorStorage{dt}([dt[] for i in 1:N_vecs])
end

function add_element!(mvs::MultiVectorStorage{T}, val::T, i::Integer) where T
    push!(mvs.data[i], val)
    return mvs
end

function Base.length(mvs::MultiVectorStorage)
    return sum(length.(mvs.data))
end

#Concats vectors into one normal vector
function Base.convert(::Type{Vector{T}}, mvs::MultiVectorStorage{T}) where T
    return reduce(vcat, mvs.data)
end