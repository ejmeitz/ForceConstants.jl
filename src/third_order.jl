export third_order_IFC

"""
Calculates analytical third order force constants for a pair potential (e.g. Lennard Jones).
For pair potentials the only non-zero terms will be self-terms (e.g. i,i,i) and terms where
i = j or j = k (e.g. 1,1,2).
"""
function third_order_IFC(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    @assert all(pot.r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    N_atoms = n_atoms(sys)
    Ψ = zeros(D*N_atoms,D*N_atoms,D*N_atoms)

    unique_indices = with_replacement_combinations((1:N_atoms), 3)
    is_not_self_term = (i,j,k) -> !(i == j && i == k)
    not_all_unique = (i,j,k) -> (i == j || i == k || j == k)

    r = zeros(D)*unit(sys.atoms.position[1][1])
    for idx in unique_indices #hard to parllelize with this iterator
        i,j,k = idx
        if is_not_self_term(i,j,k) && not_all_unique(i,j,k)
            # println(i," ",j," ",k)

            if i == j
                r .= sys.atoms.position[i] .- sys.atoms.position[k]
            elseif i == k
                r .= sys.atoms.position[i] .- sys.atoms.position[j]
            elseif j == k
                r .= sys.atoms.position[j] .- sys.atoms.position[i]
            end
            
            r .= nearest_mirror(r, sys.box_sizes_SC)
            dist = norm(r)

            if dist < pot.r_cut
                for α in range(1,D)
                    for β in range(1,D)
                        for γ in range(1,D)
                            ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(k-1) + γ
    
                            val = -ustrip(ϕ₃(pot, dist, r, α, β, γ))
                            Ψ[ii,jj,kk] = val; Ψ[ii,kk,jj] = val
                            Ψ[jj,kk,ii] = val; Ψ[jj,ii,kk] = val
                            Ψ[kk,ii,jj] = val; Ψ[kk,jj,ii] = val
                        end
                    end
                end
            end
        end
    end

    #Acoustic Sum Rule (technically re-calculating the whole loop above here)
    Threads.@threads :dynamic for i in range(1,N_atoms)
        for α in range(1,D)
            for β in range(1,D)
                for γ in range(1,D)
                    ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ
                    Ψ[ii,jj,kk] = ϕ₃_self(sys, pot, i, α, β, γ)
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
δₖ(x,y) = ==(x,y)

function ϕ₃(pot::PairPotential, r_norm, rᵢⱼ, α, β, γ)
    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    Φ′′′ = potential_third_deriv(pot, r_norm)

    return (rᵢⱼ[α]*rᵢⱼ[β]*rᵢⱼ[γ]/(r_norm^3)) * (Φ′′′ - (3*Φ′′/r_norm) + (3*Φ′/(r_norm^2))) +
        ((rᵢⱼ[α]*δₖ(β,γ) + rᵢⱼ[β]*δₖ(γ,α) + rᵢⱼ[γ]*δₖ(α,β))/(r_norm^2))*(Φ′′ - (Φ′/r_norm))

end


#Just to vectorized version
function ϕ₃_self(sys::SuperCellSystem, pot::PairPotential, i, α, β, γ)
    N_atoms = n_atoms(sys)
    value = 0.0
    
    #Loop all atoms in system except atom i
    for j in range(1,N_atoms)
        if j != i
            r_ij = position(sys, i) .- position(sys, j)
            r_ij = nearest_mirror(r_ij, sys.box_sizes_SC)
            dist = norm(r_ij)
            if dist < pot.r_cut
                value += ustrip(ϕ₃(pot, dist, r_ij, α, β, γ))
            end
        end
    end

    return value

end