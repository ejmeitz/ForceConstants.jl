export third_order, save_third_order

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

    #Acoustic Sum Rule
    Threads.@threads :dynamic for i in range(1,N_atoms)
        for α in range(1,D)
            for β in range(1,D)
                for γ in range(1,D)
                    ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ
                    Ψ[ii,jj,kk] = -1*sum(Ψ[ii, β:D:end, γ:D:end])
                    # Ψ[ii,jj,kk] = -1*sum(Ψ[α:D:end, jj, γ:D:end])
                    # Ψ[ii,jj,kk] = -1*sum(Ψ[α:D:end, β:D:end, kk])
                        Ψ[ii,jj,kk] = -1*sum(Ψ[ii, jj, γ:D:end])
                end
            end
        end
    end

    #Acoustic Sum Rule -- slow version, just debugging
    # Threads.@threads :dynamic for i in range(1,N_atoms)
    #     for α in range(1,D)
    #         for β in range(1,D)
    #             for γ in range(1,D)
    #                 ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ
    #                 Ψ[ii,jj,kk] = ϕ₃_self(sys, pot, i, α, β, γ)
    #             end
    #         end
    #     end
    # end



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

function save_third_order(Ψ, N_atoms, D, outpath; filename = "third_order", fmt = :txt)
    @assert (N_atoms*D) == size(Ψ)[1] "Incorrect dimensions"

    if fmt == :txt
        filepath = joinpath(outpath, filename*".txt")
        f = open(filepath, "w")
        for i in range(1, N_atoms)
            for j in range(1, N_atoms)
                for k in range(1,N_atoms)
                    for α in range(1,D)
                        for β in range(1,D)
                            for γ in range(1,D)
                                ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(k-1) + γ
                                if ustrip(Ψ[ii,jj,kk]) != 0.0
                                    println(f, "$i $α $j $β $k $γ $(ustrip(Ψ[ii,jj,kk]))")
                                end
                            end
                        end
                    end
                end
            end
        end
        close(f)
    elseif fmt == :HDF5
        throw(ArgumentError("Not implemented yet, $(fmt)"))
    elseif fmt == :JLD2
        throw(ArgumentError("Not implemented yet, $(fmt)"))
    else
        throw(ArgumentError("Invalid format, $(fmt)"))
    end
end

#Just to vectorized version
# function ϕ₃_self(sys::SuperCellSystem, pot::PairPotential, i, α, β, γ)
#     N_atoms = n_atoms(sys)
#     value = 0.0
    
#     #Loop all atoms in system except atom i
#     for j in range(1,N_atoms)
#         if j != i
#             r_ij = position(sys, i) .- position(sys, j)
#             r_ij = nearest_mirror(r_ij, sys.box_sizes_SC)
#             dist = norm(r_ij)
#             if dist < pot.r_cut
#                 value += ustrip(ϕ₃(pot, dist, r_ij, α, β, γ))
#             end
#         end
#     end

#     return value

# end