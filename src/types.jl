export SuperCellSystem, Potential, PairPotential, ThreeBodyPotential,
    masses, mass, position, positions, charges, charge, n_atoms, n_atoms_per_uc,
    SparseForceConstants, DenseForceConstants, F3_val, F4_val

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
abstract type FC_val{V,I} end

struct DenseForceConstants{O,V,U,T} <: AbstractForceConstants{O}
    values::Array{V,O}
    units::U
    tol::T
end

Base.size(fc::DenseForceConstants) = size(fc.values)
Base.getindex(fc::DenseForceConstants{O}, idxs::Vararg{Integer, O}) where O = fc.values[idxs...]
Base.getindex(fc::DenseForceConstants{O}, idxs::CartesianIndex{O}) where O = fc.values[idxs]

n_modes(fc::DenseForceConstants) = size(fc)[1]

#type aliases
const SecondOrderForceConstants{V,U,T} = DenseForceConstants{2,V,U,T}
const ThirdOrderForceConstants{V,U,T} = DenseForceConstants{3,V,U,T}
const FourthOrderForceConstants{V,U,T} = DenseForceConstants{4,V,U,T}

struct F4_val{V,I} <: FC_val{V,I} #Int16 is usually fine and helps with SIMD
    val::V
    i::I
    j::I
    k::I
    l::I
end

struct F3_val{V,I} <: FC_val{V,I} #Int16 is usually fine and helps with SIMD
    val::V
    i::I
    j::I
    k::I
end

value(fcv::FC_val) = fcv.val
idx(fc::F3_val) = [fc.i,fc.j,fc.k]
idx(fc::F4_val) = [fc.i,fc.j,fc.k,fc.l]

function FC_val(val, idxs::Vararg{Integer, 3})
    return F3_val{typeof(val), eltype(idxs)}(val, idxs...)
end

function FC_val(val, idxs::Vararg{Integer, 4})
    return F4_val{typeof(val), eltype(idxs)}(val, idxs...)
end


struct SparseForceConstants{O,V,U,T} <: AbstractForceConstants{O}
    values::StructArray{<:FC_val{V,I}} #1 entry per DoF
    units::U
    tol::T
end

Base.length(fc::SparseForceConstants) = sum(length.(fc.values))

#type aliases
const SparseSecondOrder{V,U,T} = SparseForceConstants{2,V,U,T}
const SparseThirdOrder{V,U,T} = SparseForceConstants{3,V,U,T}
const SparseFourthOrder{V,U,T} = SparseForceConstants{4,V,U,T}

function SparseForceConstants(vals::Vector{FC_val{V, O}}, units, tol) where {V,O}
    return SparseForceConstants{O, V, typeof(units), typeof(tol)}(StructArray(vals), units, tol)
end

# #Construct SparseForceConstants object from dense ForceConstants object 
function SparseForceConstants(fc::DenseForceConstants{O,V,U,T}; nthreads::Integer = Threads.nthreads()) where {O,V,U,T}
   
    num_nonzero = sum(fc.values .!= 0.0)
    sa = StructArray{FC_val{V,O}}(undef, (num_nonzero,))

    count = 1
    for idx in CartesianIndices(fc.values)
        if fc[idx] != 0.0
            sa[count] = FC_val(fc[idx], Tuple(idx)...)
            count += 1
        end
    end

    return SparseForceConstants{O,V,U,T}(sa, fc.units, fc.tol)
end

#&NEED TO TEST THIS
# function sort_by_mode(sfc::SparseForceConstants{O,V}, n_dof::Integer) where {O,V}
#     data = Vector{StructArray{FC_val{V,O}}}(undef, n_dof)

#     for fc in sfc.values
#         data[fc.idxs[1]] = fc
#     end

#     return data
# end

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

