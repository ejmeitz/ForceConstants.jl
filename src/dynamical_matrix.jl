export dynamicalMatrix

function dynamicalMatrix(sys::UnitCellSystem{D}, pot::PairPotential, k_point::Float64) where D
    return nothing
end

#Long term add support for non-pair potentials
function dynamicalMatrix(sys::SuperCellSystem{D}, pot::PairPotential) where D
    N_atoms = n_atoms(sys)
    dynmat = zeros(D*N_atoms,D*N_atoms)

    for i in range(1,N_atoms)
        for j in range(i+1,N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.box_sizes)
            dist = norm(rᵢⱼ)

            for α in range(1,D)
                for β in range(1,D)

                    #Calculate dynmat index
                    ii = D*(i-1) + α
                    jj = D*(j-1) + β

                    if dist < pot.r_cut
                        dynmat[ii,jj] = -ustrip(ϕ₂(pot, dist, rᵢⱼ, α, β)) #in 1D rᵢⱼ is the same as the magnitude of rᵢⱼ
                        dynmat[jj,ii] = dynmat[i,j]
                    else
                        dynmat[ii,jj] = 0.0
                        dynmat[jj,ii] = 0.0
                    end

                end
             end
        end
    end

    #Acoustic Sum Rule -- Loop D x D block matricies on diagonal of dynmat
    for i in range(1, N_atoms)
        for α in range(1,D)
            for β in range(1,D)
                ii = D*(i-1) + α
                jj = D*(i-1) + β # i == j because we're on diagonal
                dynmat[ii,jj] = -1*sum(dynmat[ii,β:D:end])
            end
        end
    end

    #Mass Weight
    for i in range(1,N_atoms)
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

    #Add final units to dynamical matrix
    dynmat *= unit(pot.ϵ / pot.σ^2 / mass(sys,1))

    return dynmat
end

function nearest_mirror(r_ij,L)
    r_x = r_ij[1]; r_y = r_ij[2]; r_z = r_ij[3]
        
    if r_x > L/2
        r_x = r_x - L
    elseif r_x < -L/2
        r_x = r_x + L
    end
        
    if r_y > L/2
        r_y = r_y - L
    elseif r_y < -L/2
        r_y = r_y + L  
    end
        
    if r_z > L/2
        r_z = r_z - L
    elseif r_z < -L/2
        r_z = r_z + L
    end
    
    return [r_x,r_y,r_z] 
end


function ϕ₂(pot, r_norm, rᵢⱼ, α, β)
    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    return (α == β) ? ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm))) + (Φ′/r_norm) :
                    ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm)))
end


function nearest_mirror(r_ij::Float64, box_sizes::AbstractVector{Int})
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