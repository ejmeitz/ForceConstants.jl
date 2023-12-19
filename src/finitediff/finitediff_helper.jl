
function energy_loop(pot::PairPotential, posns, box_sizes, N_atoms, r_cut)

    U_total = 0.0*energy_unit(pot)

    for i in range(1,N_atoms)
        for j in range(i+1, N_atoms)             
            r = posns[i] .- posns[j]
            nearest_mirror!(r, box_sizes)
            dist = norm(r)
            if dist < r_cut
                U_total += potential(pot,dist)
            end
        end
    end

    return U_total

end

#Expects positions as Vector{Vector}
function energy_loop(pot::StillingerWeberSilicon, posns, box_sizes, N_atoms, r_cut)

    U_total = 0.0*energy_unit(pot)
    # forces = zeros(N_atoms, 3)*energy_unit(pot)/length_unit(pot)

    rᵢⱼ = similar(posns[1])
    rᵢₖ = similar(posns[1])
    r_cut_sq = pot.r_cut * pot.r_cut #avoid doing sqrts
    for i in range(1,N_atoms)
        for j in range(1,N_atoms)

            if i != j
                rᵢⱼ .= posns[i] .- posns[j]
                nearest_mirror!(rᵢⱼ, box_sizes)
                dist_ij_sq = sum(x -> x^2, rᵢⱼ)

                if dist_ij_sq < r_cut_sq

                    if i < j
                        U_total += pair_potential(pot, sqrt(dist_ij_sq))
                        # forces[i,:] += 0
                        # forces[j,:] -= 0
                    end

                    for k in range(j+1, N_atoms)
                        if i != k
                            rᵢₖ .= posns[i] .- posns[k]
                            nearest_mirror!(rᵢₖ, box_sizes)
                            dist_ik_sq = sum(x -> x^2, rᵢₖ)
                            
                            if dist_ik_sq < r_cut_sq
                                U_total += three_body_potential(pot, rᵢⱼ, rᵢₖ, sqrt(dist_ij_sq), sqrt(dist_ik_sq))
                            end
                        end
                    end
                end
            end

        end
        
    end
    return U_total
end

#* slower than ^^
function energy_loop2(pot::StillingerWeberSilicon, posns, box_sizes, N_atoms, r_cut)

    U_total = 0.0*energy_unit(pot)

    rᵢⱼ = similar(posns[1])
    rᵢₖ = similar(posns[1])
    rⱼₖ = similar(posns[1])
    r_cut_sq = pot.r_cut * pot.r_cut #avoid doing sqrts
    for i in range(1, N_atoms)
        for j in range(i+1, N_atoms)

            rᵢⱼ .= posns[i] .- posns[j]
            nearest_mirror!(rᵢⱼ, box_sizes)
            dist_ij_sq = sum(x -> x^2, rᵢⱼ)

            if dist_ij_sq < r_cut_sq
                U_total += pair_potential(pot, sqrt(dist_ij_sq))
            end
            
            for k in range(j+1, N_atoms)
                rᵢₖ .= posns[i] .- posns[k]
                nearest_mirror!(rᵢₖ, box_sizes)
                dist_ik_sq = sum(x -> x^2, rᵢₖ)

                rⱼₖ .= rᵢₖ .- rᵢⱼ
                # nearest_mirror!(rⱼₖ, box_sizes)
                dist_jk_sq = sum(x -> x^2, rⱼₖ)
                
                if dist_ij_sq < r_cut_sq && dist_ik_sq < r_cut_sq
                    dist_ij = sqrt(dist_ij_sq); dist_ik = sqrt(dist_ik_sq)
                    cosθⱼᵢₖ = dot(rᵢⱼ, rᵢₖ) / (dist_ij * dist_ik)
                    U_total += h_SW(pot, dist_ij, dist_ik, cosθⱼᵢₖ)
                end
                if dist_ij_sq < r_cut_sq && dist_jk_sq < r_cut_sq
                    dist_ij = sqrt(dist_ij_sq); dist_jk = sqrt(dist_jk_sq)
                    cosθᵢⱼₖ = dot(-rᵢⱼ, rⱼₖ) / (dist_ij * dist_jk)
                    U_total += h_SW(pot, dist_ij, dist_jk, cosθᵢⱼₖ)
                end
                if dist_ik_sq < r_cut_sq && dist_jk_sq < r_cut_sq
                    dist_jk = sqrt(dist_jk_sq); dist_ik = sqrt(dist_ik_sq)
                    cosθᵢₖⱼ = dot(rⱼₖ, rᵢₖ) / (dist_jk * dist_ik)
                    U_total += h_SW(pot, dist_ik, dist_jk, cosθᵢₖⱼ)
                end
                
            end
        end
    end
    return U_total
end



# Calculates force on atom j in β direction
function force_loop_j(pot::PairPotential, posns, force_unit, r_cut, box_sizes, N_atoms, j, β)

    D = length(box_sizes)
    Fⱼᵦ = zeros(D)*force_unit

    for i in range(1,N_atoms)
        if i != j            
            r = posns[i] .- posns[j]
            nearest_mirror!(r, box_sizes)
            dist = norm(r)
            if dist < r_cut
                F_mag = force(pot, dist)
                r_hat = r / dist 
                F_ij = F_mag*r_hat

                #Update forces and potential
                Fⱼᵦ += F_ij[β]
            end
        end
    end

    return Fⱼᵦ

end
