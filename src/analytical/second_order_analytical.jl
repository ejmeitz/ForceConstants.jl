export second_order_analytical

function second_order_analytical(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)

    #Loop block matricies above diagonal (not including diagonal)
    for i in range(1,N_atoms)
        for j in range(i+1,N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
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

    IFC2 = apply_tols!(IFC2,tol)

    return DenseForceConstants(IFC2, energy_unit(pot) / length_unit(pot)^2, tol)

end

function ϕ₂(pot::PairPotential, r_norm, r_jk_j′k′, α, β)
    rᵢⱼ = r_jk_j′k′

    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    return (α == β) ? ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm))) + (Φ′/r_norm) :
                    ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm)))
end