export fourth_order_IFC

"""
Calculates analytical fourth order force constants for a pair potential (e.g. Lennard Jones).
For pair potentials the only non-zero terms will be self-terms (e.g. i,i,i,i) and terms where
i = j or j = k (e.g. 1,1,2) and terms where i = j = k (e.g. (i,i,i,j)). All terms that have
fully unique indicies will be 0.
"""
function fourth_order_IFC(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    @assert all(pot.r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    N_atoms = n_atoms(sys)
    χ = zeros(D*N_atoms,D*N_atoms,D*N_atoms,D*N_atoms)

    unique_indices = with_replacement_combinations((1:N_atoms), 4)
    is_not_self_term = (i,j,k,l) -> !(i == j && i == k && i == l)
    not_all_unique = (i,j,k,l) -> (length([i,j,k,l]) != length(unique([i,j,k,l])))

    r = zeros(D)*unit(sys.atoms.position[1][1])
    for idx in unique_indices #hard to parllelize with this iterator
        i,j,k,l = idx
        if is_not_self_term(i,j,k,l) && not_all_unique(i,j,k)

            # REDO
            # if i == j
            #     r .= sys.atoms.position[i] .- sys.atoms.position[k]
            # elseif i == k
            #     r .= sys.atoms.position[i] .- sys.atoms.position[j]
            # elseif j == k
            #     r .= sys.atoms.position[j] .- sys.atoms.position[i]
            # end
            
            r .= nearest_mirror(r, sys.box_sizes_SC)
            dist = norm(r)

            if dist < pot.r_cut
                for α in range(1,D)
                    for β in range(1,D)
                        for γ in range(1,D)
                            for δ in range(1,D)
                                ii = D*(i-1) + α; jj = D*(j-1) + β
                                kk = D*(k-1) + γ; ll = D*(l-1) + δ
        
                                val = -ustrip(ϕ₄(pot, dist, r, α, β, γ, δ))
                                for idx2 in permutations([ii,jj,kk,ll],4)
                                    χ[idx2...] = val
                                end
                            end
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
                    for δ in range(1,D)
                        ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ; ll = D*(i-1) + δ
                        χ[ii,jj,kk,ll] = ϕ₄_self(sys, pot, i, α, β, γ, δ)
                    end
                end
            end
        end
    end

    #Apply tolerances
    χ[abs.(χ) .< tol] .= 0.0

    #Give proper units
    χ *= unit(pot.ϵ / pot.σ^4)

    return χ
end

function fourth_order_sparse(sys::SuperCellSystem{D}, pot::PairPotential) where D

end

#Kronicker Delta
δₖ(x,y) = ==(x,y)

function ϕ₄(pot::PairPotential, r_norm, rᵢⱼ, α, β, γ, δ)
    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    Φ′′′ = potential_third_deriv(pot, r_norm)
    Φ′′′′ = potential_fourth_deriv(pot, r_norm)

    v = sum([rᵢⱼ[α]*rᵢⱼ[β], rᵢⱼ[α]*rᵢⱼ[γ], rᵢⱼ[α]*rᵢⱼ[δ], rᵢⱼ[β]*rᵢⱼ[γ], rᵢⱼ[β]*rᵢⱼ[δ], rᵢⱼ[γ]*rᵢⱼ[δ]] .*
         [δₖ(γ,δ),δₖ(β,δ),δₖ(β,γ),δₖ(α,δ),δₖ(α,γ),δₖ(α,β)])

    return (rᵢⱼ[α]*rᵢⱼ[β]*rᵢⱼ[γ]*rᵢⱼ[δ]/(r_norm^4)) * (Φ′′′′ - (6*Φ′′′/r_norm) + (15*Φ′′/(r_norm^2)) - (15*Φ′/(r_norm^3))) +
        (v/(r_norm^3))*(Φ′′′ - (3*Φ′′/r_norm) + (3*Φ′/(r_norm^2))) + 
        (δₖ(α, β)*δₖ(γ,δ) + δₖ(α,γ)*δₖ(β,δ) + δₖ(α,δ)*δₖ(β,γ))*((Φ′′-(Φ′/r))/(r_norm^2))

end


#Just to vectorized version
function ϕ₄_self(sys::SuperCellSystem, pot::PairPotential, i, α, β, γ, δ)
    N_atoms = n_atoms(sys)
    value = 0.0
    
    #Loop all atoms in system except atom i
    for j in range(1,N_atoms)
        if j != i
            r_ij = position(sys, i) .- position(sys, j)
            r_ij = nearest_mirror(r_ij, sys.box_sizes_SC)
            dist = norm(r_ij)
            if dist < pot.r_cut
                value += ustrip(ϕ₃(pot, dist, r_ij, α, β, γ, δ))
            end
        end
    end

    return value

end