export second_order, second_order!

function second_order(sys::SuperCellSystem{D}, pot::PairPotential, 
    calc::AnalyticalCalculator) where D

    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)
    r_cut_sq = calc.r_cut*calc.r_cut

    #Loop block matricies above diagonal (not including diagonal)
    for i in range(1,N_atoms)
        for j in range(i+1,N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
            dist_sq = sum(x -> x^2, rᵢⱼ)

            if dist_sq < r_cut_sq
                dist = norm(rᵢⱼ)
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

    IFC2 = apply_tols!(IFC2, calc.tol)

    return DenseForceConstants(IFC2, energy_unit(pot) / length_unit(pot)^2, calc.tol)

end


function second_order!(IFC2, sys::SuperCellSystem{D}, pot::PairPotential,
    calc::AnalyticalCalculator) where D

    N_atoms = n_atoms(sys)
    @assert size(IFC2) == (D*N_atoms,D*N_atoms)

    #Loop block matricies above diagonal (not including diagonal)
    for i in range(1,N_atoms)
        for j in range(i+1,N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
            dist_sq = sum(x -> x^2, rᵢⱼ)

            if dist_sq < r_cut_sq
                dist = norm(rᵢⱼ)
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

    return IFC2

end

function ϕ₂(pot::PairPotential, r_norm, r_jk_j′k′, α, β)
    rᵢⱼ = r_jk_j′k′

    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    return (α == β) ? ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm))) + (Φ′/r_norm) :
                    ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm)))
end