export fourth_order_sparse

"""
Calculates analytical fourth order force constants for a pair potential (e.g. Lennard Jones).
For pair potentials the only non-zero terms will be self-terms (e.g. i,i,i,i) and terms where
i = j and k = l (e.g. (i,i,j,j)) and terms where i = j = k (e.g. (i,i,i,j)). All terms that have
more than 2 unique indices will be 0.
"""
function fourth_order_sparse(sys::SuperCellSystem{D}, pot::PairPotential,
         calc::AnalyticalCalculator; nthreads::Integer = Threads.nthreads()) where D

    @assert all(r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    
    N_atoms = n_atoms(sys)

    if (N_atoms < nthreads)
        @warn "More threads than atoms, using N_atoms as number of threads"
        nthreads = N_atoms
    end

    if (N_atoms <= typemax(Int16)) #* this probably causes type instabilitiy
        int_type = Int16
    else
        int_type = Int32
    end

    r_cut_sq = calc.r_cut*calc.r_cut
    
    #Create storage for each thread
    χ_arrays = MultiVectorStorage(F4_val{Float64,int_type}, nthreads)

    atoms_per_thd = cld(N_atoms,nthreads) #* DONT LIKE THIS, BAD BALANCING, USE ITERATORS.PARTITION


    Base.@sync for thd in 1:nthreads
        Base.Threads.@spawn begin
            χ_thd = χ_arrays.data[thd] #Get local storage for force constants
            start_idx = int_type((thd-1)*atoms_per_thd + 1)
            end_idx = int_type(min(thd*atoms_per_thd, N_atoms))

            r = zeros(D)*unit(sys.atoms.position[1][1]) #pre-allocate
            for i in range(start_idx, end_idx)
                for j in range(int_type(i+1),N_atoms)

                    # i,i,i,j terms
                    r .= sys.atoms.position[i] .- sys.atoms.position[j]
                    nearest_mirror!(r, sys.box_sizes_SC)
                    dist_sq = sum(x -> x^2, r)

                    if dist_sq < r_cut_sq
                        dist = norm(r)
                        for α in range(1,D)
                            for β in range(1,D)
                                for γ in range(1,D)
                                    for δ in range(1,D)
                                        val = ustrip(ϕ₄(pot, dist, r, α, β, γ, δ))
                                        
                                        if abs(val) > calc.tol
                                            # i,i,j,j terms
                                            ii = int_type(D*(i-1) + α); jj = int_type(D*(i-1) + β)
                                            kk = int_type(D*(j-1) + γ); ll = int_type(D*(j-1) + δ)
                                            set_terms_fourth_order!(χ_thd, ii, jj, kk, ll, val)

                                            # i,i,i,j terms
                                            ii = int_type(D*(i-1) + α); jj = int_type(D*(i-1) + β)
                                            kk = int_type(D*(i-1) + γ); ll = int_type(D*(j-1) + δ)
                                            set_terms_fourth_order!(χ_thd, ii, jj, kk, ll, -val)

                                        end
                                        
                                    end
                                end
                            end
                        end
                    end

                    #i,j,j,j
                    r .= sys.atoms.position[j] .- sys.atoms.position[i]
                    nearest_mirror!(r, sys.box_sizes_SC)
                    dist_sq = sum(x -> x^2, r)

                    if dist_sq < r_cut_sq
                        dist = norm(r)
                        for α in range(1,D)
                            for β in range(1,D)
                                for γ in range(1,D)
                                    for δ in range(1,D)
                                        val = -ustrip(ϕ₄(pot, dist, r, α, β, γ, δ))
                                        
                                        if abs(val) > calc.tol
                                            # j,j,j,i terms
                                            ii = int_type(D*(i-1) + α); jj = int_type(D*(j-1) + β)
                                            kk = int_type(D*(j-1) + γ); ll = int_type(D*(j-1) + δ)
                                            set_terms_fourth_order!(χ_thd, ii, jj, kk, ll, val)

                                            #&  j,j,i,i terms just same as above
                                        end
                                        
                                    end
                                end
                            end
                        end 
                    end

                                                    
                end
            end
        end
    end

    #Concat arrays calculated by separate threads
    χ = reduce(vcat, χ_arrays.data)

    #Acoustic Sum Rule (i,i,i,i terms)
    #* Accumulate these in the loop above
    Threads.@threads for i in range(1,N_atoms)
        for α in range(1,D)
            for β in range(1,D)
                for γ in range(1,D)
                    for δ in range(1,D)
                        ii = int_type(D*(i-1) + α); jj = int_type(D*(i-1) + β);
                        kk = int_type(D*(i-1) + γ); ll = int_type(D*(i-1) + δ)
                        val = ϕ₄_self(sys, pot, i, α, β, γ, δ, r_cut)
                        if abs(val) > calc.tol
                            push!(χ,FC_val(val, ii, jj, kk, ll))
                        end
                    end
                end
            end
        end
    end

    #Give proper units
    χ_unit = energy_unit(pot) / length_unit(pot)^4

    return SparseForceConstants(χ, χ_unit, calc.tol)
end

function set_terms_fourth_order!(χ, ii, jj, kk, ll, val)
    #Terms starting with ii
    push!(χ, FC_val(val,ii,jj,kk,ll))
    push!(χ, FC_val(val,ii,kk,jj,ll))
    push!(χ, FC_val(val,ii,jj,ll,kk))
    push!(χ, FC_val(val,ii,jj,kk,ll))
    push!(χ, FC_val(val,ii,ll,jj,kk))
    push!(χ, FC_val(val,ii,kk,ll,jj))

    #Terms starting with jj
    push!(χ, FC_val(val,jj,ii,kk,ll))
    push!(χ, FC_val(val,jj,kk,ii,ll))
    push!(χ, FC_val(val,jj,ii,ll,kk))
    push!(χ, FC_val(val,jj,ll,ii,kk))
    push!(χ, FC_val(val,jj,ll,kk,ii))
    push!(χ, FC_val(val,jj,kk,ll,ii))

    #Terms starting with kk
    push!(χ, FC_val(val,kk,ii,jj,ll))
    push!(χ, FC_val(val,kk,jj,ii,ll))
    push!(χ, FC_val(val,kk,ii,ll,jj))
    push!(χ, FC_val(val,kk,ll,ii,jj))
    push!(χ, FC_val(val,kk,ll,jj,ii))
    push!(χ, FC_val(val,kk,jj,ll,ii))

    #Terms starting with ll
    push!(χ, FC_val(val,ll,ii,jj,kk))
    push!(χ, FC_val(val,ll,jj,ii,kk))
    push!(χ, FC_val(val,ll,ii,kk,jj))
    push!(χ, FC_val(val,ll,kk,ii,jj))
    push!(χ, FC_val(val,ll,kk,jj,ii))
    push!(χ, FC_val(val,ll,jj,kk,ii))

    return χ
end

function potential_derivs4(pot::PairPotential, r)
    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    Φ′′′ = potential_third_deriv(pot, r_norm)
    Φ′′′′ = potential_fourth_deriv(pot, r_norm)

    return Φ′, Φ′′,Φ′′′,Φ′′′′
end

function ϕ₄(pot::PairPotential, r_norm, rᵢⱼ, α, β, γ, δ)
    Φ′ = potential_first_deriv(pot, r_norm) #* pull out these are independent of abcd
    Φ′′ = potential_second_deriv(pot, r_norm)
    Φ′′′ = potential_third_deriv(pot, r_norm)
    Φ′′′′ = potential_fourth_deriv(pot, r_norm)

    v = rᵢⱼ[α]*rᵢⱼ[β]*δₖ(γ,δ) + rᵢⱼ[α]*rᵢⱼ[γ]*δₖ(β,δ) + rᵢⱼ[α]*rᵢⱼ[δ]*δₖ(β,γ) + rᵢⱼ[β]*rᵢⱼ[γ]*δₖ(α,δ) +
         rᵢⱼ[β]*rᵢⱼ[δ]*δₖ(α,γ) + rᵢⱼ[γ]*rᵢⱼ[δ]*δₖ(α,β)


    return ((rᵢⱼ[α]*rᵢⱼ[β]*rᵢⱼ[γ]*rᵢⱼ[δ]/(r_norm^4)) * (Φ′′′′ - (6*Φ′′′/r_norm) + (15*Φ′′/(r_norm^2)) - (15*Φ′/(r_norm^3)))) +
        ((v/(r_norm^3))*(Φ′′′ - (3*Φ′′/r_norm) + (3*Φ′/(r_norm^2)))) + 
        (((δₖ(α, β)*δₖ(γ,δ)) + (δₖ(α,γ)*δₖ(β,δ)) + (δₖ(α,δ)*δₖ(β,γ)))*((Φ′′-(Φ′/r_norm))/(r_norm^2)))

end


function ϕ₄_self(sys::SuperCellSystem, pot::PairPotential, i, α, β, γ, δ, r_cut)
    N_atoms = n_atoms(sys)
    value = 0.0
    
    #Loop all atoms in system except atom i
    for j in range(1,N_atoms)
        if j != i
            r_ij = position(sys, i) .- position(sys, j)
            r_ij = nearest_mirror(r_ij, sys.box_sizes_SC)
            dist = norm(r_ij)
            if dist < r_cut
                value += ustrip(ϕ₄(pot, dist, r_ij, α, β, γ, δ))
            end
        end
    end

    return value

end


function mass_weight_fourth_order!(ifc4::SparseForceConstants{4}, masses)
    mass_unit = unit(masses[1])
    masses = ustrip.(masses)

    for (i,v) in enumerate(eachindex(ifc4.values))
        ifc4.values[i].val /= sqrt.(masses[v.idxs])
    end

    ifc4.units /= mass_unit
    return ifc4
end


