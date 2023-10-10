export fourth_order_sparse, fourth_order_finite_diff

"""
Calculates analytical fourth order force constants for a pair potential (e.g. Lennard Jones).
For pair potentials the only non-zero terms will be self-terms (e.g. i,i,i,i) and terms where
i = j and k = l (e.g. (i,i,j,j)) and terms where i = j = k (e.g. (i,i,i,j)). All terms that have
more than 2 unique indices will be 0.
"""
function fourth_order_sparse(sys::SuperCellSystem{D}, pot::PairPotential,
         tol; float_type = Float32, r_cut = pot.r_cut, nthreads::Integer = Threads.nthreads()) where D

    @assert all(r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    
    N_atoms = n_atoms(sys)

    if (N_atoms < nthreads)
        @warn "More threads than atoms, using N_atoms as number of threads"
        nthreads = N_atoms
    end

    #Create storage for each thread
    χ_arrays = MultiVectorStorage(FC_val{float_type,4}, nthreads)

    atoms_per_thd = cld(N_atoms,nthreads) #* DONT LIKE THIS, BAD BALANCING

    Base.@sync for thd in 1:nthreads
        Base.Threads.@spawn begin
            χ_thd = χ_arrays.data[thd] #Get local storage for force constants
            start_idx = (thd-1)*atoms_per_thd + 1
            end_idx = min(thd*atoms_per_thd, N_atoms)

            r = zeros(D)*unit(sys.atoms.position[1][1]) #pre-allocate
            for i in range(start_idx, end_idx)
                for j in range(i+1,N_atoms)

                    # i,i,i,j terms
                    r .= sys.atoms.position[i] .- sys.atoms.position[j]
                    nearest_mirror!(r, sys.box_sizes_SC)
                    dist = norm(r)

                    if dist < r_cut
                        for α in range(1,D)
                            for β in range(1,D)
                                for γ in range(1,D)
                                    for δ in range(1,D)
                                        val = float_type(ustrip(ϕ₄(pot, dist, r, α, β, γ, δ)))
                                        
                                        if abs(val) > tol
                                            # i,i,j,j terms
                                            ii = D*(i-1) + α; jj = D*(i-1) + β
                                            kk = D*(j-1) + γ; ll = D*(j-1) + δ
                                            set_terms_fourth_order!(χ_thd, ii, jj, kk, ll, val)

                                            # i,i,i,j terms
                                            ii = D*(i-1) + α; jj = D*(i-1) + β
                                            kk = D*(i-1) + γ; ll = D*(j-1) + δ
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

                    if dist < r_cut
                        for α in range(1,D)
                            for β in range(1,D)
                                for γ in range(1,D)
                                    for δ in range(1,D)
                                        val = -float_type(ustrip(ϕ₄(pot, dist, r, α, β, γ, δ)))
                                        
                                        if abs(val) > tol
                                            # j,j,j,i terms
                                            ii = D*(i-1) + α; jj = D*(j-1) + β
                                            kk = D*(j-1) + γ; ll = D*(j-1) + δ
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
    χ = convert(Vector{FC_val{float_type, 4}}, χ_arrays)

    #Acoustic Sum Rule (i,i,i,i terms) #& how to parallelize this part?
    for i in range(1,N_atoms)
        for α in range(1,D)
            for β in range(1,D)
                for γ in range(1,D)
                    for δ in range(1,D)
                        ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ; ll = D*(i-1) + δ
                        val = float_type(ϕ₄_self(sys, pot, i, α, β, γ, δ, r_cut))
                        if abs(val) > tol
                            push!(χ,FC_val(val, ii, jj, kk, ll))
                        end
                    end
                end
            end
        end
    end

    #Give proper units
    χ_unit = unit(pot.ϵ / pot.σ^4)

    return χ
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

#https://arxiv.org/pdf/2104.04895.pdf
#Eqn 9
function fourth_order_finite_diff(sys::SuperCellSystem{3}, pot::PairPotential, atom_idxs, cartesian_idxs;
     r_cut = pot.r_cut, h = 0.04*0.5291772109u"Å")

    h = uconvert(unit(pot.σ), h)

    @assert length(atom_idxs) == 4
    @assert length(cartesian_idxs) == 4
    @assert !all(atom_idxs .== atom_idxs[1]) "Cannot check self terms with finite difference"

    N_atoms = n_atoms(sys)

    combos = [[h,h,h],[h,h,-h],[h,-h,h],[h,-h,-h],[-h,h,h],[-h,h,-h],[-h,-h,h],[-h,-h,-h]]
    
    force_unit = zero(force(pot,1u"Å")) 
    forces = zeros(8)*force_unit

    #Make mutable #& change SimpleCrystals to not use SVector
    posns = [Vector(a) for a in positions(sys)]

    for (c,combo) in enumerate(combos)

        posns[atom_idxs[1]][1] += combo[1]
        posns[atom_idxs[2]][2] += combo[2]
        posns[atom_idxs[3]][3] += combo[3]

        l = atom_idxs[4]
        force_val = zero(force(pot,1u"Å"))  
        #Calculate force on atom l
        for i in range(1,N_atoms)
            if i != l               
                r = posns[i] .- posns[l]
                nearest_mirror!(r, sys.box_sizes_SC)
                dist = norm(r)
                
                #Make sure mirrored particle is in cuttoff
                if dist < r_cut
                    #Get force, potential with modified potential/ force function
                    F_mag = force(pot, dist)
                    
                    r_hat = r / dist 
                    F_ij = F_mag.*r_hat

                    #Update forces and potential
                    force_val += F_ij[cartesian_idxs[4]] #& IS THSI SIGN RIGHT??
                end
            end
        end

        forces[c] = force_val

        #Un-modify
        posns[atom_idxs[1]][1] -= combo[1]
        posns[atom_idxs[2]][2] -= combo[2]
        posns[atom_idxs[3]][3] -= combo[3]
    end

    return (1/(8*(h^3)))*(forces[1] - forces[2] - forces[3] + forces[4] - forces[5] + forces[6] + forces[7] - forces[8])
end