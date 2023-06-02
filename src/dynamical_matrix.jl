export dynamicalMatrix, second_order_IFC

### Unit Cell ###
function dynamicalMatrix(sys::UnitCellSystem{D}, pot::PairPotential, k_point::SVector{D}, tol) where D
    @assert all(pot.r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    
    atoms_per_unit_cell = n_atoms_per_uc(sys)
    dynmat = zeros(ComplexF64, D*atoms_per_unit_cell, D*atoms_per_unit_cell)

    dynamicalMatrix_UC_Helper!(sys, pot, dynmat, k_point)

    #Remove entries smaller than tol
    dynmat[abs.(dynmat) .< tol] .= 0.0

    #Add final units to dynamical matrix
    dynmat *= unit(pot.ϵ / pot.σ^2 / mass(sys,1))

    return dynmat
end

function dynamicalMatrix_UC_Helper!(sys::UnitCellSystem{2}, pot::PairPotential, dynmat, k_point)
    atoms_per_unit_cell = n_atoms_per_uc(sys)

    #Loop block matricies on and above diagonal
    for i in range(1, atoms_per_unit_cell)
        # Position of atom i in base unitcell
        r_i0 = position(sys, (1,1), i)
        for j in range(i, atoms_per_unit_cell)
            for α in range(1,2)
                for β in range(1,2)
                    #Calculate dynmat index
                    ii = 2*(i-1) + α
                    jj = 2*(j-1) + β

                    #Sum over all unit-cells
                    for uc_idx in keys(sys.atoms)
                        #Position of atom j in unitcell k
                        r_jk = position(sys, uc_idx, j)

                        r_i0_jk = r_i0 .- r_jk
                        r_i0_jk = nearest_mirror(r_i0_jk, sys.box_sizes_SC)
                        dist = norm(r_i0_jk)

                        #Self terms
                        if j == i && uc_idx == (1,1) #dist is 0 so no exp term
                            dynmat[ii,jj] += ϕ₂_self(sys, pot, α, β, r_i0, i)
                        elseif dist < pot.r_cut
                            ϕ₂_val = -ustrip(ϕ₂(pot, dist, r_i0_jk, α, β))
                            exp_part = exp(-im*dot(ustrip(k_point), ustrip(r_i0_jk)))
                            if imag(exp_part) < 1e-7
                                dynmat[ii,jj] += ϕ₂_val*real(exp_part)
                            else
                                raise(error("Imaginary part too large: $(exp_part), power: $(-dot(ustrip(k_point), ustrip(r_i0_jk)))"))
                            end
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

                        r_i0_jk = r_i0 .- r_jk
                        r_i0_jk = nearest_mirror(r_i0_jk, sys.box_sizes_SC)
                        dist = norm(r_i0_jk)

                        #Self terms
                        if j == i && uc_idx == (1,1,1) #dist is 0 so no exp term
                            dynmat[ii,jj] += ϕ₂_self(sys, pot, α, β, r_i0, i)
                        elseif dist < pot.r_cut
                             dynmat[ii,jj] += -ustrip(ϕ₂(pot, dist, r_i0_jk, α, β))*
                                exp(-im*dot(ustrip(k_point), ustrip(r_i0_jk)))
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
    @assert all(pot.r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    N_atoms = n_atoms(sys)

    #reuse storage from IFC2 calculation
    dynmat = second_order_IFC(sys, pot, tol, false)

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

function second_order_IFC(sys::SuperCellSystem{D}, pot::PairPotential, tol, add_units) where D
    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)

    #Loop block matricies above diagonal (not including diagonal)
    for i in range(1,N_atoms)
        for j in range(i+1,N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.box_sizes_SC)
            dist = norm(rᵢⱼ)

            if dist < pot.r_cut
                for α in range(1,D)
                    for β in range(1,D)

                        #Calculate IFC2 index
                        ii = D*(i-1) + α
                        jj = D*(j-1) + β

                        IFC2[ii,jj] = -ustrip(ϕ₂(pot, dist, rᵢⱼ, α, β))
                        IFC2[jj,ii] = IFC2[ii,jj]
                    end
                end
            end
        end
    end

    #Acoustic Sum Rule -- Loop D x D block matricies on diagonal of IFC2
    for i in range(1, N_atoms) # index of block matrix
        for α in range(1,D)
            for β in range(1,D)
                ii = D*(i-1) + α
                jj = D*(i-1) + β # i == j because we're on diagonal
                IFC2[ii,jj] = -1*sum(IFC2[ii, β:D:end])
            end
        end
    end

    #Remove entries smaller than tol
    IFC2[abs.(IFC2) .< tol] .= 0.0

    #Add final units to dynamical matrix
    if add_units
        IFC2 *= unit(pot.ϵ / pot.σ^2)
    end

    return IFC2

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
function ϕ₂_self(sys::UnitCellSystem{D}, pot::PairPotential, α, β, r_i0, i) where D
    atoms_per_unit_cell = n_atoms_per_uc(sys)
    value = 0.0
    
    #Loop all atoms in system
    for uc_idx in keys(sys.atoms)
        for j in range(1, atoms_per_unit_cell)
            if !(uc_idx == tuple(ones(D)...) && i == j) #skip atom i0
                r_jk = position(sys, uc_idx, j)

                r_i0_jk = r_i0 .- r_jk
                r_i0_jk = nearest_mirror(r_i0_jk, sys.box_sizes_SC)
                dist = norm(r_i0_jk)
                
                if dist < pot.r_cut
                    value += ustrip(ϕ₂(pot, dist, r_i0_jk, α, β))
                end
            end
        end
    end

    return value
end

#could remove allocations here and just make it a ! function
function nearest_mirror(r_ij, box_sizes)
    r_x = r_ij[1]; r_y = r_ij[2]; r_z = r_ij[3]
    L_x, L_y, L_z = box_sizes
    if r_x > L_x/2
        r_x -= L_x
    elseif r_x < -L_x/2
        r_x += L_x
    end
        
    if r_y > L_y/2
        r_y -= L_y
    elseif r_y < -L_y/2
        r_y += L_y  
    end
        
    if r_z > L_z/2
        r_z -= L_z
    elseif r_z < -L_z/2
        r_z += L_z
    end
    
    return [r_x,r_y,r_z] 
end