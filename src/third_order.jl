export third_order_IFC, third_order_test, mass_weight_sparsify_third_order, mass_weight_third_order!, F3_val

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
    for idx in unique_indices #TODO PARALLELIZE THIS, HOW TO STILL USE ITERATOR? just for i <= j <= k
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
    Threads.@threads for i in range(1,N_atoms)
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
    Ψ = apply_tols!(Ψ, tol)

    #Give proper units
    Ψ_unit = unit(pot.ϵ / pot.σ^3)

    return ForceConstants(Ψ, Ψ_unit, tol)
end

function mass_weight_third_order!(Ψ::ThirdOrderForceConstants, masses::AbstractVector)
    N_modes = size(Ψ)[1]
    D = Int(N_modes/length(masses))
    N_atoms = length(masses)

    @assert D ∈ [1,2,3] 

    mass_unit = unit(masses[1])
    masses = ustrip.(masses)

    for i in 1:N_atoms
        for j in 1:N_atoms
            for k in 1:N_atoms
                for α in 1:D
                    for β in 1:D
                        for γ in 1:D
                            ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(k-1) + γ
                            if Ψ[ii,jj,kk] != 0
                                Ψ.values[ii,jj,kk] /= sqrt(masses[i]*masses[j]*masses[k])
                            end
                        end
                    end
                end
            end
        end
    end

    Ψ.units /= mass_unit
    return Ψ
end

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


"""
Mass weights the force constant matrix such that element i,j,k is divided by
sqrt(m_i * m_j * m_k). This is useful when converting to modal coupling constants.
The force constants are also returned in sparse format as a vector of values and indices.
"""
function mass_weight_sparsify_third_order(Ψ::ThirdOrderForceConstants, masses::AbstractVector)

    N_modes = size(Ψ)[1]
    D = Int(N_modes/length(masses))
    N_atoms = length(masses)

    @assert D ∈ [1,2,3] 

    mass_unit = unit(masses[1])
    masses = ustrip.(masses)

    num_nonzero = sum(Ψ.values .!= 0.0)
    Ψ_non_zero_mw = Vector{F3_val}(undef,(num_nonzero,))

    count = 1
    for i in 1:N_atoms
        for j in 1:N_atoms
            for k in 1:N_atoms
                for α in 1:D
                    for β in 1:D
                        for γ in 1:D
                            ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(k-1) + γ
                            if Ψ[ii,jj,kk] != 0
                                Ψ_non_zero_mw[count] = F3_val(ii, jj, kk, Ψ[ii,jj,kk]/sqrt(masses[i]*masses[j]*masses[k]))
                                count += 1
                            end
                        end
                    end
                end
            end
        end
    end
    @assert (count - 1) == num_nonzero


    return SparseForceConstants(Ψ_non_zero_mw, Ψ.units/mass_unit, Ψ.tol)
end


function third_order_test(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    @assert all(pot.r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    N_atoms = n_atoms(sys)
    Ψ = zeros(D*N_atoms, D*N_atoms, D*N_atoms)


    # pot is pair potential only loop atomic pairs
    Threads.@threads for i in range(1,N_atoms)
        r = zeros(D)*unit(sys.atoms.position[1][1]) #pre-allocate per thread
        for j in range(i+1,N_atoms)

            #i,i,j terms
            r .= sys.atoms.position[i] .- sys.atoms.position[j]
            nearest_mirror!(r, sys.box_sizes_SC)
            dist = norm(r)

            
            if dist < pot.r_cut
                for α in range(1,D)
                    for β in range(1,D)
                        for γ in range(1,D)
                            val = -ustrip(ϕ₃(pot, dist, r, α, β, γ)) #* re-calculating derivatives ununecessarily in here
                            ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(j-1) + γ
                            Ψ[ii,jj,kk] = val; Ψ[ii,kk,jj] = val
                            Ψ[jj,kk,ii] = val; Ψ[jj,ii,kk] = val
                            Ψ[kk,ii,jj] = val; Ψ[kk,jj,ii] = val 
                        end
                    end
                end            

            end

            # #j,j,i term
            r .= sys.atoms.position[j] .- sys.atoms.position[i]
            nearest_mirror!(r, sys.box_sizes_SC)
            dist = norm(r)

            if dist < pot.r_cut
                for α in range(1,D)
                    for β in range(1,D)
                        for γ in range(1,D)
                            val = -ustrip(ϕ₃(pot, dist, r, α, β, γ)) 
                            ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(j-1) + γ
                            Ψ[ii,jj,kk] = val; Ψ[ii,kk,jj] = val
                            Ψ[jj,kk,ii] = val; Ψ[jj,ii,kk] = val
                            Ψ[kk,ii,jj] = val; Ψ[kk,jj,ii] = val 
                        end
                    end
                end            

            end
        end
    end

    Ψ = ASR!(Ψ, N_atoms, D)

    #Apply tolerances
    Ψ = apply_tols!(Ψ, tol)

    #Give proper units
    Ψ_unit = unit(pot.ϵ / pot.σ^3)

    return DenseForceConstants(Ψ, Ψ_unit, tol)
    
end


function ASR!(ifc3::Array{T,3}, N_atoms, D) where T
    #Loop all self terms
    Threads.@threads for i in range(1,N_atoms)
        for α in range(1,D)
            for β in range(1,D)
                for γ in range(1,D)
                    ii_self = D*(i-1) + α; jj_self = D*(i-1) + β; kk_self = D*(i-1) + γ 

                    # Loop atoms k with this α, β, γ
                    for k in range(1,N_atoms)
                        if k != i
                            ii = D*(i-1) + α; jj = D*(ji-1) + β; kk = D*(k-1) + γ
                            ifc3[ii_self,jj_self,kk_self] -= ifc3[ii,jj,kk] 
                        end
                    end

                end
            end
        end
    end

    return ifc3
end