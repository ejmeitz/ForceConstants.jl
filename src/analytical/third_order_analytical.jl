export third_order, third_order!, mass_weight_third_order!

function third_order(sys::SuperCellSystem{D}, pot::PairPotential,
    calc::AnalyticalCalculator) where D

    @assert all(pot.r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    N_atoms = n_atoms(sys)
    Ψ = zeros(D*N_atoms, D*N_atoms, D*N_atoms)

    r_cut_sq = calc.r_cut*calc.r_cut

    # pot is pair potential only loop atomic pairs
    Threads.@threads for i in range(1,N_atoms)
        r = zeros(D)*unit(sys.atoms.position[1][1]) #pre-allocate per thread
        for j in range(i+1,N_atoms)

            #i,i,j terms
            r .= sys.atoms.position[i] .- sys.atoms.position[j]
            nearest_mirror!(r, sys.box_sizes_SC)
            dist_sq = sum(x -> x^2, r)

            
            if dist_sq < r_cut_sq
                dist = norm(r)
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
            dist_sq = sum(x -> x^2, r)

            if dist_sq < r_cut_sq
                dist = norm(r)
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

    #Self terms
    Ψ = ASR!(Ψ, N_atoms, D)

    #Apply tolerances
    Ψ = apply_tols!(Ψ, calc.tol)

    #Give proper units
    Ψ_unit = energy_unit(pot) / length_unit(pot)^3

    return DenseForceConstants(Ψ, Ψ_unit, calc.tol)
    
end

function third_order!(Ψ, sys::SuperCellSystem{D}, pot::PairPotential,
    calc::AnalyticalCalculator) where D

    @assert all(pot.r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    N_atoms = n_atoms(sys)
    @assert size(Ψ) == (D*N_atoms, D*N_atoms, D*N_atoms)

    r_cut_sq = calc.r_cut*calc.r_cut

    # pot is pair potential only loop atomic pairs
    Threads.@threads for i in range(1,N_atoms)
        r = zeros(D)*unit(sys.atoms.position[1][1]) #pre-allocate per thread
        for j in range(i+1,N_atoms)

            #i,i,j terms
            r .= sys.atoms.position[i] .- sys.atoms.position[j]
            nearest_mirror!(r, sys.box_sizes_SC)
            dist_sq = sum(x -> x^2, r)

            
            if dist_sq < r_cut_sq
                dist = norm(r)
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
            dist_sq = sum(x -> x^2, r)
            
            #ASR
            if dist_sq < r_cut_sq
                dist = norm(r)
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

    #Self terms
    Ψ = ASR!(Ψ, N_atoms, D)

    return Ψ
    
end

#* untested
# function third_order_sparse(sys::SuperCellSystem{D}, pot::PairPotential,
#     tol; r_cut = pot.r_cut, nthreads::Integer = Threads.nthreads()) where D

#     throw(error("Not yet implemented"))

#     @assert all(r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    
#     N_atoms = n_atoms(sys)

#     if (N_atoms < nthreads)
#         @warn "More threads than atoms, using N_atoms as number of threads"
#         nthreads = N_atoms
#     end

#     #Can help with SIMD later
#     if (N_atoms <= typemax(Int16))
#         int_type = Int16
#     else
#         int_type = Int32
#     end

#     #Create storage for each thread
#     Ψ_arrays = MultiVectorStorage(F3_val{Float64,int_type}, nthreads)
#     atoms_per_thd = cld(N_atoms,nthreads) #* DONT LIKE THIS, BAD BALANCING, USE ITERATORS.PARTITION

#     Base.@sync for thd in 1:nthreads
#         Base.Threads.@spawn begin
#             Ψ_thd = Ψ_arrays.data[thd] #Get local storage for force constants
#             start_idx = (thd-1)*atoms_per_thd + 1
#             end_idx = min(thd*atoms_per_thd, N_atoms)

#             r = zeros(D)*unit(sys.atoms.position[1][1]) #pre-allocate
#             for i in range(start_idx, end_idx)
#                 for j in range(i+1,N_atoms)

#                     r .= sys.atoms.position[i] .- sys.atoms.position[j]
#                     nearest_mirror!(r, sys.box_sizes_SC)
#                     dist = norm(r)
        
                    
#                     if dist < pot.r_cut
#                         for α in range(1,D)
#                             for β in range(1,D)
#                                 for γ in range(1,D)
#                                     val = -ustrip(ϕ₃(pot, dist, r, α, β, γ)) #* re-calculating derivatives ununecessarily in here
#                                     ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(j-1) + γ
#                                     push!(Ψ_thd, FC_val(val,ii,jj,kk))
#                                     push!(Ψ_thd, FC_val(val,ii,kk,jj))
#                                     push!(Ψ_thd, FC_val(val,jj,kk,ii))
#                                     push!(Ψ_thd, FC_val(val,jj,ii,kk))
#                                     push!(Ψ_thd, FC_val(val,kk,ii,jj))
#                                     push!(Ψ_thd, FC_val(val,kk,jj,ii))
#                                 end
#                             end
#                         end            
        
#                     end
        
#                     # #j,j,i term
#                     r .= sys.atoms.position[j] .- sys.atoms.position[i]
#                     nearest_mirror!(r, sys.box_sizes_SC)
#                     dist = norm(r)
        
#                     if dist < pot.r_cut
#                         for α in range(1,D)
#                             for β in range(1,D)
#                                 for γ in range(1,D)
#                                     val = -ustrip(ϕ₃(pot, dist, r, α, β, γ)) 
#                                     ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(j-1) + γ
#                                     push!(Ψ_thd, FC_val(val,ii,jj,kk))
#                                     push!(Ψ_thd, FC_val(val,ii,kk,jj))
#                                     push!(Ψ_thd, FC_val(val,jj,kk,ii))
#                                     push!(Ψ_thd, FC_val(val,jj,ii,kk))
#                                     push!(Ψ_thd, FC_val(val,kk,ii,jj))
#                                     push!(Ψ_thd, FC_val(val,kk,jj,ii))
#                                 end
#                             end
#                         end            
        
#                     end
#                 end
#             end
#         end
#     end

#     Ψ = convert(Vector{FC_val{Float64, 3}}, Ψ_arrays)

#     #Acoustic Sum Rule (i,i,i,i terms) #& how to parallelize this part?
#     #* Can keep running track of these terms when generating others
#     for i in range(1,N_atoms)
#         for α in range(1,D)
#             for β in range(1,D)
#                 for γ in range(1,D)
#                     ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ
#                     val = ϕ₃_self(sys, pot, i, α, β, γ)
#                     if abs(val) > tol
#                         push!(Ψ ,FC_val(val, ii, jj, kk))
#                     end
#                 end
#             end
#         end
#     end


#     #Give proper units
#     Ψ_unit = energy_unit(pot) / length_unit(pot)^3

#     return SparseForceConstants(Ψ, Ψ_unit, tol)

# end

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
