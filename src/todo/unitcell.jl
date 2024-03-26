# This code should work but was more written as a way for me to understand wahts going on
# It does not have any optimizations and should be just as expensive as the supercell calculation.

# struct UnitCellSystem{D,L} <: System
#     atoms::Dict{Tuple,StructArray{Atom}}
#     num_unit_cells::SVector{D,Int}
#     box_sizes_UC::L
#     box_sizes_SC::L
# end

# function UnitCellSystem(crystal::Crystal{D}, num_unit_cells::SVector{D,Int}) where D
#     @assert sum(crystal.N_unit_cells) == D "Crystal must be a single unit cell"

#     unit_cell_indices = Iterators.product([1:num_unit_cells[i] for i in range(1,D)]...)

#     positions = SimpleCrystals.position(crystal)
#     masses = SimpleCrystals.atomic_mass(crystal)
#     charges = getindex.(crystal.atoms, :charge)
#     box_sizes = norm.(bounding_box(crystal))
#     box_sizes_SC = norm.(num_unit_cells .* bounding_box(crystal))
#     #Generate data for all atoms in all unit cells
#     atoms = Dict()

#     for uc_idx in unit_cell_indices
#         positions_new = []
#         unitcell_origin = [box_sizes[1]*(uc_idx[1]-1), box_sizes[2]*(uc_idx[2]-1), box_sizes[3]*(uc_idx[3]-1)]
#         for i in range(1, length(crystal))
#             push!(positions_new, positions[i] .+ unitcell_origin)
#         end
#         atoms[uc_idx] = StructArray{Atom}((position = positions_new, mass = masses, charge = charges))
#     end

#     return UnitCellSystem{D, typeof(box_sizes)}(atoms, num_unit_cells, box_sizes, box_sizes_SC)
# end

# mass(sys::UnitCellSystem, atom_idx::Int) = sys.atoms[(1,1,1)].mass[atom_idx]
# position(sys::UnitCellSystem, uc_idx::Tuple, atom_idx::Int) = sys.atoms[uc_idx].position[atom_idx]
# charge(sys::UnitCellSystem, atom_idx::Int) = sys.atoms[(1,1,1)].charge[atom_idx]
# n_atoms(sys::UnitCellSystem) = length(sys.atoms)*length(sys.atoms[(1,1,1)])
# n_atoms_per_uc(sys::UnitCellSystem) = length(sys.atoms[(1,1,1)])


# ### Unit Cell ###
# function dynamicalMatrix(sys::UnitCellSystem{D}, pot::PairPotential, k_point::SVector{D}, tol) where D
#     @assert all(pot.r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    
#     atoms_per_unit_cell = n_atoms_per_uc(sys)
#     dynmat = zeros(ComplexF64, D*atoms_per_unit_cell, D*atoms_per_unit_cell)

#     dynamicalMatrix_UC_Helper!(sys, pot, dynmat, k_point)

#     #Remove entries smaller than tol
#     for i in eachindex(dynmat)
#         if abs(dynmat[i]) < tol
#             dynmat[i] = 0.0
#         end
#     end

#     #Add final units to dynamical matrix
#     dynmat_unit = unit(pot.ϵ / pot.σ^2 / mass(sys,1))

#     return SecondOrderMatrix(dynmat, dynmat_unit, tol)
# end

# function dynamicalMatrix_UC_Helper!(sys::UnitCellSystem{2}, pot::PairPotential, dynmat, k_point)
#     atoms_per_unit_cell = n_atoms_per_uc(sys)

#     #Loop block matricies on and above diagonal
#     for i in range(1, atoms_per_unit_cell)
#         # Position of atom i in base unitcell
#         r_i0 = position(sys, (1,1), i)
#         for j in range(i, atoms_per_unit_cell)
#             for α in range(1,2)
#                 for β in range(1,2)
#                     #Calculate dynmat index
#                     ii = 2*(i-1) + α
#                     jj = 2*(j-1) + β

#                     #Sum over all unit-cells
#                     for uc_idx in keys(sys.atoms)
#                         #Position of atom j in unitcell k
#                         r_jk = position(sys, uc_idx, j)

#                         r_i0_jk = r_i0 .- r_jk
#                         r_i0_jk = nearest_mirror(r_i0_jk, sys.box_sizes_SC)
#                         dist = norm(r_i0_jk)

#                         #Self terms
#                         if j == i && uc_idx == (1,1) #dist is 0 so no exp term
#                             dynmat[ii,jj] += ϕ₂_self(sys, pot, α, β, r_i0, i)
#                         elseif dist < pot.r_cut
#                             ϕ₂_val = -ustrip(ϕ₂(pot, dist, r_i0_jk, α, β))
#                             exp_part = exp(-im*dot(ustrip(k_point), ustrip(r_i0_jk)))
#                             if imag(exp_part) < 1e-7
#                                 dynmat[ii,jj] += ϕ₂_val*real(exp_part)
#                             else
#                                 raise(error("Imaginary part too large: $(exp_part), power: $(-dot(ustrip(k_point), ustrip(r_i0_jk)))"))
#                             end
#                         end

#                     end
#                     #Enforce symmetry & Mass Weighting
#                     dynmat[ii,jj] /= ustrip(sqrt(mass(sys,i)*mass(sys,j)))
#                     dynmat[jj,ii] = dynmat[ii,jj]                    
#                 end
#             end
#         end
#     end

#     return dynmat
# end


# function dynamicalMatrix_UC_Helper!(sys::UnitCellSystem{3}, pot::PairPotential, dynmat, k_point)
#     atoms_per_unit_cell = n_atoms_per_uc(sys)

#     #Loop block matricies on and above diagonal
#     for i in range(1, atoms_per_unit_cell)
#         # Position of atom i in base unitcell
#         r_i0 = position(sys, (1,1,1), i)
#         for j in range(i, atoms_per_unit_cell)

#             for α in range(1,3)
#                 for β in range(1,3)
#                     #Calculate dynmat index
#                     ii = 3*(i-1) + α
#                     jj = 3*(j-1) + β

#                     #Sum over all unit-cells
#                     for uc_idx in keys(sys.atoms)
#                         #Position of atom j in unitcell k
#                         r_jk = position(sys, uc_idx, j)

#                         r_i0_jk = r_i0 .- r_jk
#                         r_i0_jk = nearest_mirror(r_i0_jk, sys.box_sizes_SC)
#                         dist = norm(r_i0_jk)

#                         #Self terms
#                         if j == i && uc_idx == (1,1,1) #dist is 0 so no exp term
#                             dynmat[ii,jj] += ϕ₂_self(sys, pot, α, β, r_i0, i)
#                         elseif dist < pot.r_cut
#                              dynmat[ii,jj] += -ustrip(ϕ₂(pot, dist, r_i0_jk, α, β))*
#                                 exp(-im*dot(ustrip.(k_point), ustrip.(r_i0_jk)))
#                         end

#                     end
#                     #Enforce symmetry & Mass Weighting
#                     dynmat[ii,jj] /= ustrip(sqrt(mass(sys,i)*mass(sys,j)))
#                     dynmat[jj,ii] = dynmat[ii,jj]                    
#                 end
#             end
#         end
#     end

#     return dynmat
# end


#Not efficient but it works -- reuse the ϕ₂ terms???
# function ϕ₂_self(sys::UnitCellSystem{D}, pot::PairPotential, α, β, r_i0, i) where D
#     atoms_per_unit_cell = n_atoms_per_uc(sys)
#     value = 0.0
    
#     #Loop all atoms in system
#     for uc_idx in keys(sys.atoms)
#         for j in range(1, atoms_per_unit_cell)
#             if !(uc_idx == tuple(ones(D)...) && i == j) #skip atom i0
#                 r_jk = position(sys, uc_idx, j)

#                 r_i0_jk = r_i0 .- r_jk
#                 r_i0_jk = nearest_mirror(r_i0_jk, sys.box_sizes_SC)
#                 dist = norm(r_i0_jk)
                
#                 if dist < pot.r_cut
#                     value += ustrip(ϕ₂(pot, dist, r_i0_jk, α, β))
#                 end
#             end
#         end
#     end

#     return value
# end
