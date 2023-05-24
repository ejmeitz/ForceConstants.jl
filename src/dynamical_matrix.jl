function DynamicalMatrix(sys::UnitCellSystem{D}, pot::PairPotential, k_point::Float64) where {D <: Integer}

end

#Long term add support for non-pair potentials
function DynamicalMatrix(sys::SuperCellSystem{D}, pot::PairPotential) where {D <: Integer}
    N_atoms = n_atoms(sys)
    dynmat = zeros(D*N_atoms,D*N_atoms)

    for i in range(1,N_atoms)
        for j in range(i+1,N_atoms)

            rᵢⱼ = norm(sys.atoms.position[i] .- sys.atoms.position[j])
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.L)
            dist = norm(rᵢⱼ)

            for α in range(1,3)
                for β in range(1,3)

                    #Calculate dynmat index
                    ii = 3*(i-1) + α
                    jj = 3*(j-1) + β

                    if dist < r_cut
                        dynmat[ii,jj] = -ϕ₂(pot,dist, rᵢⱼ, α, β) #in 1D rᵢⱼ is the same as the magnitude of rᵢⱼ
                        dynmat[jj,ii] = dynmat[i,j]
                    else
                        dynmat[ii,jj] = 0.0
                        dynmat[jj,ii] = 0.0
                    end

                end
             end
        end
    end

    #Acoustic Sum Rule
    for i in range(1, D*N_atoms)
        dynmat[i,i] = -1*sum(dynmat[i,:])
    end

    #Mass Weight
    for i in range(1,D*N_atoms)
        for j in range(1, D*N_atoms)
            dynmat[i,j] *=  (1/sqrt(mass(sys,i)*mass(sys,j)))
        end
    end
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