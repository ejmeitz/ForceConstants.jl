export fourth_order_sparse

"""
Calculates analytical fourth order force constants for a pair potential (e.g. Lennard Jones).
For pair potentials the only non-zero terms will be self-terms (e.g. i,i,i,i) and terms where
i = j and k = l (e.g. (i,i,j,j)) and terms where i = j = k (e.g. (i,i,i,j)). All terms that have
fully unique indicies will be 0.
"""
function fourth_order_sparse(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    @assert all(pot.r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    N_atoms = n_atoms(sys)
    Χ = FC_val[]

    for i in range(1,N_atoms)
        r = zeros(D)*unit(sys.atoms.position[1][1])

        for j in range(i+1,N_atoms)

            # i,i,i,j terms
            r .= sys.atoms.position[i] .- sys.atoms.position[j]
            nearest_mirror!(r, sys.box_sizes_SC)
            dist = norm(r)

            if dist < pot.r_cut
                for α in range(1,D)
                    for β in range(1,D)
                        for γ in range(1,D)
                            for δ in range(1,D)
                                val = ustrip(ϕ₄(pot, dist, r, α, β, γ, δ))
                                
                                if val != 0.0
                                    # i,i,j,j terms
                                    ii = D*(i-1) + α; jj = D*(i-1) + β
                                    kk = D*(j-1) + γ; ll = D*(j-1) + δ
                                    set_terms_fourth_order!(χ, ii, jj, kk, ll, val)

                                    # i,i,i,j terms
                                    ii = D*(i-1) + α; jj = D*(i-1) + β
                                    kk = D*(i-1) + γ; ll = D*(j-1) + δ
                                    set_terms_fourth_order!(χ, ii, jj, kk, ll, -val)

                                end
                                
                            end
                        end
                    end
                end
            end

            #i,j,j,j
            r .= sys.atoms.position[j] .- sys.atoms.position[i]
            nearest_mirror!(r, sys.box_sizes_SC)

            if dist < pot.r_cut
                for α in range(1,D)
                    for β in range(1,D)
                        for γ in range(1,D)
                            for δ in range(1,D)
                                val = -ustrip(ϕ₄(pot, dist, r, α, β, γ, δ)) 
                                
                                if val != 0.0
                                     # j,j,j,i terms
                                     ii = D*(i-1) + α; jj = D*(j-1) + β
                                     kk = D*(j-1) + γ; ll = D*(j-1) + δ
                                     set_terms_fourth_order!(χ, ii, jj, kk, ll, val)

                                     #&  j,j,i,i terms just same as above
                                end
                                
                            end
                        end
                    end
                end
            end

                                               
        end
    end

    #Acoustic Sum Rule (i,i,i,i terms)
    Threads.@threads for i in range(1,N_atoms)
        for α in range(1,D)
            for β in range(1,D)
                for γ in range(1,D)
                    for δ in range(1,D)
                        ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ; ll = D*(i-1) + δ
                        append!(χ,FC_val(ϕ₄_self(sys, pot, i, α, β, γ, δ), ii, jj, kk, ll))
                    end
                end
            end
        end
    end

    #Apply tolerances #& IS THIS REALLY EXPENSIVE WITH DELETIONS?
    χ = apply_tols!(χ,tol)

    #Give proper units
    χ *= unit(pot.ϵ / pot.σ^4)

    return χ
end

function set_terms_fourth_order!(χ::AbstractVector{FC_val}, ii, jj, kk, ll, val)
    #Terms starting with ii
    append!(χ, FC_val(val,ii,jj,kk,ll))
    append!(χ, FC_val(val,ii,kk,jj,ll))
    append!(χ, FC_val(val,ii,jj,ll,kk))
    append!(χ, FC_val(val,ii,jj,kk,ll))
    append!(χ, FC_val(val,ii,ll,jj,kk))
    append!(χ, FC_val(val,ii,kk,ll,jj))

    #Terms starting with jj
    append!(χ, FC_val(val,jj,ii,kk,ll))
    append!(χ, FC_val(val,jj,kk,ii,ll))
    append!(χ, FC_val(val,jj,ii,ll,kk))
    append!(χ, FC_val(val,jj,ll,ii,kk))
    append!(χ, FC_val(val,jj,ll,kk,ii))
    append!(χ, FC_val(val,jj,kk,ll,ii))

    #Terms starting with kk
    append!(χ, FC_val(val,kk,ii,jj,ll))
    append!(χ, FC_val(val,kk,jj,ii,ll))
    append!(χ, FC_val(val,kk,ii,ll,jj))
    append!(χ, FC_val(val,kk,ll,ii,jj))
    append!(χ, FC_val(val,kk,ll,jj,ii))
    append!(χ, FC_val(val,kk,jj,ll,ii))

    #Terms starting with ll
    append!(χ, FC_val(val,ll,ii,jj,kk))
    append!(χ, FC_val(val,ll,jj,ii,kk))
    append!(χ, FC_val(val,ll,ii,kk,jj))
    append!(χ, FC_val(val,ll,kk,ii,jj))
    append!(χ, FC_val(val,ll,kk,jj,ii))
    append!(χ, FC_val(val,ll,jj,kk,ii))

    return χ
end


function ϕ₄(pot::PairPotential, r_norm, rᵢⱼ, α, β, γ, δ)
    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    Φ′′′ = potential_third_deriv(pot, r_norm)
    Φ′′′′ = potential_fourth_deriv(pot, r_norm)

    v = rᵢⱼ[α]*rᵢⱼ[β]*δₖ(γ,δ) + rᵢⱼ[α]*rᵢⱼ[γ]*δₖ(β,δ) + rᵢⱼ[α]*rᵢⱼ[δ]*δₖ(β,γ) + rᵢⱼ[β]*rᵢⱼ[γ]*δₖ(α,δ) +
         rᵢⱼ[β]*rᵢⱼ[δ]*δₖ(α,γ) + rᵢⱼ[γ]*rᵢⱼ[δ]*δₖ(α,β)


    return ((rᵢⱼ[α]*rᵢⱼ[β]*rᵢⱼ[γ]*rᵢⱼ[δ]/(r_norm^4)) * (Φ′′′′ - (6*Φ′′′/r_norm) + (15*Φ′′/(r_norm^2)) - (15*Φ′/(r_norm^3)))) +
        ((v/(r_norm^3))*(Φ′′′ - (3*Φ′′/r_norm) + (3*Φ′/(r_norm^2)))) + 
        (((δₖ(α, β)*δₖ(γ,δ)) + (δₖ(α,γ)*δₖ(β,δ)) + (δₖ(α,δ)*δₖ(β,γ)))*((Φ′′-(Φ′/r_norm))/(r_norm^2)))

end


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
                value += ustrip(ϕ₄(pot, dist, r_ij, α, β, γ, δ))
            end
        end
    end

    return value

end