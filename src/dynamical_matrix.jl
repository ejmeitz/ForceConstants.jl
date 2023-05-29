export dynamicalMatrix

### Unit Cell ###
function dynamicalMatrix(sys::UnitCellSystem{D}, pot::PairPotential, k_point::SVector{D}, tol) where D
    atoms_per_unit_cell = n_atoms_per_uc(sys)
    dynmat = zeros(D*atoms_per_unit_cell, D*atoms_per_unit_cell)

    dynamicalMatrix_UC_Helper!(sys, pot, dynmat, k_point)

    #Remove entries smaller than tol
    dynmat[abs.(dynmat) .< tol] .= 0.0

    #Add final units to dynamical matrix
    dynmat *= unit(pot.ϵ / pot.σ^2 / mass(sys,1))

    return dynmat
end

function dynamicalMatrix_UC_Helper!(sys::UnitCellSystem{2}, pot::PairPotential, dynmat, k_point)
    atoms_per_unit_cell = n_atoms(sys)

    #Loop block matricies above diagonal (not including diagonal)
    for i in range(1, atoms_per_unit_cell)
        # Position of atom i in base unitcell
        r_i0 = sys.atoms.position[i]
        for j in range(1, atoms_per_unit_cell)
            for α in range(1,2)
                for β in range(1,2)
                    #Calculate dynmat index
                    ii = 2*(i-1) + α
                    jj = 2*(j-1) + β

                    #Sum over all unit-cells
                    for k1 in range(1,sys.num_unit_cells[1])
                        for k2 in range(1,sys.num_unit_cells[2])

                            #Position of atom j in unitcell k
                            kth_unitcell_origin = [sys.box_sizes[1]*(k1-1), sys.box_sizes[2]*(k2-1)]
                            r_jk = kth_unitcell_origin .+ sys.atoms.position[j]

                            r_jk_i0 = r_jk .- r_i0
                            r_jk_i0 = nearest_mirror(r_jk_i0, sys.box_sizes)
                            dist = norm(r_jk_i0)

                            if dist < pot.r_cut
                                dynmat[ii,jj] += -ustrip(ϕ₂(pot, dist, r_jk_i0, α, β))*exp(im*dot(ustrip(k_point), ustrip(r_jk_i0)))
                            end
                        end
                    end

                    dynmat[ii,jj] /= ustrip(sqrt(mass(sys,i)*mass(sys,j)))
                    dynmat[jj,ii] = dynmat[ii,jj]
                end
             end
        end
    end

    return dynmat
end


function dynamicalMatrix_UC_Helper!(sys::UnitCellSystem{3}, pot::PairPotential, dynmat, k_point)
    atoms_per_unit_cell = n_atoms_per_uc(sys)

    #Loop block matricies on and above diagonal
    for i in range(1, atoms_per_unit_cell)
        # Position of atom i in base unitcell
        r_i0 = position(sys, (1,1,1), i)
        for j in range(i, atoms_per_unit_cell)
            for α in range(1,3)
                for β in range(1,3)
                    #Calculate dynmat index
                    ii = 3*(i-1) + α
                    jj = 3*(j-1) + β

                    #Sum over all unit-cells
                    for uc_idx in keys(sys.atoms)
                        #Position of atom j in unitcell k
                        r_jk = position(sys, uc_idx, j)

                        r_jk_i0 = r_jk .- r_i0
                        r_i0_jk = r_i0 .- r_jk
                        r_jk_i0 = nearest_mirror(r_jk_i0, sys.box_sizes_SC)
                        r_i0_jk = nearest_mirror(r_i0_jk, sys.box_sizes_SC)
                        dist = norm(r_jk_i0)

                        #Cross terms
                        if dist < pot.r_cut && j != i
                            ϕ₂_val = -ustrip(ϕ₂(pot, dist, r_i0_jk, α, β))
                            dynmat[ii,jj] += ϕ₂_val*exp(im*dot(ustrip(k_point), ustrip(r_jk_i0)))
                        #Self terms
                        elseif uc_idx == (1,1,1) && j == i #dist is always 0 for self terms --> exp() term is 1
                            dynmat[ii,jj] += ϕ₂_self(sys, pot, α, β, r_i0, i)
                        end
 
                    end
                    #Enforce symmetry & Mass Weighting
                    dynmat[ii,jj] /= ustrip(sqrt(mass(sys,i)*mass(sys,j)))
                    dynmat[jj,ii] = dynmat[ii,jj]                    
                end
            end
        end
    end

    return dynmat
end

### Super Cell ###

function dynamicalMatrix(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    N_atoms = n_atoms(sys)
    dynmat = zeros(D*N_atoms,D*N_atoms)

    #Loop block matricies above diagonal (not including diagonal)
    for i in range(1,N_atoms)
        for j in range(i+1,N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.box_sizes_SC)
            dist = norm(rᵢⱼ)

            for α in range(1,D)
                for β in range(1,D)

                    #Calculate dynmat index
                    ii = D*(i-1) + α
                    jj = D*(j-1) + β

                    if dist < pot.r_cut
                        dynmat[ii,jj] = -ustrip(ϕ₂(pot, dist, rᵢⱼ, α, β))
                        dynmat[jj,ii] = dynmat[ii,jj]
                    end

                end
             end
        end
    end

    #Acoustic Sum Rule -- Loop D x D block matricies on diagonal of dynmat
    for i in range(1, N_atoms) # index of block matrix
        for α in range(1,D)
            for β in range(1,D)
                ii = D*(i-1) + α
                jj = D*(i-1) + β # i == j because we're on diagonal
                dynmat[ii,jj] = -1*sum(dynmat[ii, β:D:end])
            end
        end
    end

    #Mass Weight
    for i in range(1, N_atoms)
        for j in range(1, N_atoms)
            for α in range(1,D)
                for β in range(1,D)
                    ii = D*(i-1) + α
                    jj = D*(j-1) + β
                    dynmat[ii,jj] /=  ustrip(sqrt(mass(sys,i)*mass(sys,j)))
                end
            end
        end
    end

    #Remove entries smaller than tol
    dynmat[abs.(dynmat) .< tol] .= 0.0

    #Add final units to dynamical matrix
    dynmat *= unit(pot.ϵ / pot.σ^2 / mass(sys,1))

    return dynmat
end

### Utility Functions ###

function ϕ₂(pot::PairPotential, r_norm, r_jk_j′k′, α, β)
    rᵢⱼ = r_jk_j′k′

    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    return (α == β) ? ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm))) + (Φ′/r_norm) :
                    ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm)))
end

#Not efficient but it works -- reuse the ϕ₂ terms???
function ϕ₂_self(sys, pot, α, β, r_i0, i)
    atoms_per_unit_cell = n_atoms_per_uc(sys)
    value = 0.0
    
    #Loop all atoms in system
    for j in range(1, atoms_per_unit_cell)
        for uc_idx in keys(sys.atoms)
            if (uc_idx == (1,1,1) && i == j) #skip atom i0
                continue
            else
                r_jk = position(sys, uc_idx, j)

                r_i0_jk = r_i0 .- r_jk
                r_i0_jk = nearest_mirror(r_i0_jk, sys.box_sizes_SC)
                dist = norm(r_i0_jk)

                value += ustrip(ϕ₂(pot, dist, r_i0_jk, α, β))
            end
        end
    end

    return value
end

function nearest_mirror(r_ij, box_sizes)
    r_x = r_ij[1]; r_y = r_ij[2]; r_z = r_ij[3]
    L_x, L_y, L_z = box_sizes
    if r_x > L_x/2
        r_x = r_x - L_x
    elseif r_x < -L_x/2
        r_x = r_x + L_x
    end
        
    if r_y > L_y/2
        r_y = r_y - L_y
    elseif r_y < -L_y/2
        r_y = r_y + L_y  
    end
        
    if r_z > L_z/2
        r_z = r_z - L_z
    elseif r_z < -L_z/2
        r_z = r_z + L_z
    end
    
    return [r_x,r_y,r_z] 
end