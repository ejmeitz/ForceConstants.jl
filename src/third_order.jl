export third_order

"""
Calculates analytical third order force constants for a pair potential (e.g. Lennard Jones).
For pair potentials the only non-zero terms will be self-terms (e.g. i,i,i) and terms where
i = j or j = k (e.g. 1,1,2).
"""
function third_order(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    N_atoms = n_atoms(sys)
    Ψ = zeros(D*N_atoms,D*N_atoms,D*N_atoms)

    unique_indices = with_replacement_combinations((1:N_atoms), 3)
    is_not_self_term = (i,j,k) -> !(i == j && i == k)
    not_all_unique = (i,j,k) -> (i == j || i == k || j == k)

    for idx = unique_indices
        i,j,k = idx
        for α in range(1,D)
            for β in range(1,D)
                for γ in range(1,D)
                    if is_not_self_term(i,j,k) && not_all_unique(i,j,k)
                        ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(k-1) + γ
                        
                        if i == j
                            r = sys.atoms.position[i] .- sys.atoms.position[j]
                        elseif i == k
                            r = sys.atoms.position[i] .- sys.atoms.position[k]
                        elseif j == k
                            r = sys.atoms.position[j] .- sys.atoms.position[k]
                        end

                        r = nearest_mirror(r, sys.box_sizes_SC)
                        dist = norm(r)

                        val = -ustrip(ϕ₃(pot, dist, r, α, β, γ))
                        Ψ[ii,jj,kk] = val; Ψ[ii,kk,jj] = val
                        Ψ[jj,kk,ii] = val; Ψ[jj,ii,kk] = val
                        Ψ[kk,ii,jj] = val; Ψ[kk,jj,ii] = val
                    end
                end
            end
        end
    end

    #Acoustic Sum Rule
    for i in range(1,N_atoms)
        for α in range(1,D)
            for β in range(1,D)
                for γ in range(1,D)
                    ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ
                    Ψ[ii,jj,kk] = -1*sum(Ψ[ii, β:D:end, γ:D:end])
                end
            end
        end
    end

    #Apply tolerances
    Ψ[abs.(Ψ) .< tol] .= 0.0

    
    #Give proper units
    Ψ *= unit(pot.ϵ / pot.σ^3)

    return Ψ
end

function third_order_sparse(sys::SuperCellSystem{D}, pot::PairPotential) where D

end

#Kronicker Delta
δ(x,y) = ==(x,y)

function ϕ₃(pot::PairPotential, r_norm, rᵢⱼ, α, β, γ)
    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    Φ′′′ = potential_third_deriv(pot, r_norm)

    return (rᵢⱼ[α]*rᵢⱼ[β]*rᵢⱼ[γ]/(r_norm^3))*(Φ′′′ - (3*Φ′′/r_norm) + (3*Φ′/(r_norm^2))) +
        ((rᵢⱼ[α]*δ(β,γ) + rᵢⱼ[β]*δ(γ,α) + rᵢⱼ[γ]*δ(α,β))/(r_norm^2))*(Φ′′ - (Φ′/r_norm))

end