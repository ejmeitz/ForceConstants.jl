export second_order

function second_order(sys::SuperCellSystem{D}, pot::PairPotential,
      calc::AutoDiffCalculator) where D

    vars = make_variables(:r, D)
    r_norm = sqrt(sum(x -> x^2, vars))
    pot_symbolic = potential_nounits(pot, r_norm)
    H_symbolic = hessian(pot_symbolic, vars)
    H_exec = make_function(H_symbolic, vars)
    r_cut_sq = calc.r_cut*calc.r_cut

    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)

    for i in range(1, N_atoms)
        for j in range(i + 1, N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
            dist_ij_sq = sum(x -> x^2, rᵢⱼ)

            if dist_ij_sq < r_cut_sq
                ij_block = H_exec(ustrip.(rᵢⱼ))

                IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(j-1) + 1 : D*(j-1) + D] .= ij_block
                IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(i-1) + 1 : D*(i-1) + D] .= ij_block
            end
        end
    end

    #Acoustic Sum Rule
    Threads.@threads for i in range(1, N_atoms) # index of block matrix
        for α in range(1,D)
            for β in range(1,D)
                ii = D*(i-1) + α
                jj = D*(i-1) + β # i == j because we're on diagonal
                IFC2[ii,jj] = -1*sum(IFC2[ii, β:D:end])
            end
        end
    end

    #Taking deriv of r_ij gives you a negative sign on all terms
    # dU/dr_i_a = - dU/dr_ij_a
    IFC2 = -1 .* IFC2 

    IFC2 = apply_tols!(IFC2, calc.tol)

    return DenseForceConstants(IFC2, energy_unit(pot) / length_unit(pot)^2, calc.tol)

end

function second_order(sys::SuperCellSystem{D}, pot::StillingerWeberSilicon,
     calc::AutoDiffCalculator) where D

    @assert calc.r_cut <= pot.r_cut "For SW silicon force constant 
        cutoff must be less than potential cutoff"

    H2_exec, H3_exec_ij, H3_exec_jk, H3_exec_ik = 
        three_body_second_derivs(pot, D)
   
    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)
    r_cut_sq = calc.r_cut*calc.r_cut  

    #Loop Atomic Interactions and Add their contribution to various derivatives
    Threads.@threads for A in range(1,N_atoms)
        block = zeros(D,D)
        r_ab = similar(sys.atoms.position[1])
        rᵢₖ = similar(sys.atoms.position[1])
        nearest_j = similar(sys.atoms.position[1])
        nearest_k = similar(sys.atoms.position[1])
        r_arr = Vector{Float64}(undef, D*D)
        for B in range(1, N_atoms)
            if A != B
                #Two body term
                r_ab .= sys.atoms.position[A] .- sys.atoms.position[B]
                nearest_mirror!(r_ab, sys.box_sizes_SC)
                dist_ab_sq = sum(x -> x^2, r_ab)

                if dist_ab_sq < r_cut_sq
                    if A < B
                        block .= H2_exec(ustrip.(r_ab))
                        IFC2[D*(A-1) + 1 : D*(A-1) + D, D*(B-1) + 1 : D*(B-1) + D] .-= block
                        IFC2[D*(B-1) + 1 : D*(B-1) + D, D*(A-1) + 1 : D*(A-1) + D] .-= block
                    end

                    #Three body terms:
                    for k in range(B+1, N_atoms)
                        if A != k
                            rᵢₖ .= sys.atoms.position[A] .- sys.atoms.position[k]
                            nearest_mirror!(rᵢₖ, sys.box_sizes_SC)
                            dist_ik_sq = sum(x -> x^2, rᵢₖ)
                            
                            if dist_ik_sq < r_cut_sq
                                #ϕ₃(rᵢⱼ, rᵢₖ) contributes to ∂ri∂rj, ∂ri∂rk, ∂rj∂rk (and the symmetric versions)
                                nearest_j .= sys.atoms.position[A] .- r_ab
                                nearest_k .= sys.atoms.position[A] .- rᵢₖ
                 
                                #ij contribution
                                r_arr .= ustrip.([sys.atoms.position[A]; nearest_j; nearest_k])
                                block .= H3_exec_ij(r_arr)
                                IFC2[D*(A-1) + 1 : D*(A-1) + D, D*(B-1) + 1 : D*(B-1) + D] .+= block
                                IFC2[D*(B-1) + 1 : D*(B-1) + D, D*(A-1) + 1 : D*(A-1) + D] .+= block
                                
                                #ik contribution
                                block .= H3_exec_ik(r_arr)
                                IFC2[D*(A-1) + 1 : D*(A-1) + D, D*(k-1) + 1 : D*(k-1) + D] .+= block
                                IFC2[D*(k-1) + 1 : D*(k-1) + D, D*(A-1) + 1 : D*(A-1) + D] .+= block
                                 
                                #jk contribution
                                block .= H3_exec_jk(r_arr)
                                IFC2[D*(B-1) + 1 : D*(B-1) + D, D*(k-1) + 1 : D*(k-1) + D] .+= block
                                IFC2[D*(k-1) + 1 : D*(k-1) + D, D*(B-1) + 1 : D*(B-1) + D] .+= block
                            end
                        end
                    end
                end
            end

        end
    end

    #Acoustic Sum Rule
    Threads.@threads for i in range(1, N_atoms) # index of block matrix
        for α in range(1,D)
            for β in range(1,D)
                ii = D*(i-1) + α
                jj = D*(i-1) + β # i == j because we're on diagonal
                IFC2[ii,jj] = -1*sum(IFC2[ii, β:D:end])
            end
        end
    end

    IFC2 = apply_tols!(IFC2, calc.tol)

    return DenseForceConstants(IFC2, energy_unit(pot) / length_unit(pot)^2, calc.tol)

end


# function second_order_check(sys::SuperCellSystem{D}, pot::StillingerWeberSilicon,
#     calc::AutoDiffCalculator) where D
#     #loop through blocks in IFC2
#     for i in range(1, N_atoms)
#         for j in range(i + 1, N_atoms)

#             #Derivative of two body term is nonzero only for r_ij term
#             rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
#             rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
#             dist_ij_sq = sum(x -> x^2, rᵢⱼ)

#             if dist_ij_sq < r_cut_sq
#                 ij_block = H_exec(ustrip.(rᵢⱼ))

#                 IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(j-1) + 1 : D*(j-1) + D] .= ij_block
#                 IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(i-1) + 1 : D*(i-1) + D] .= ij_block
#             end


#             #Derivative of three body term is non-zero in multiple places

#             #When i,j = 1st,2nd atoms in triplet
#             for k in range(j+1)
#                 rᵢₖ = sys.atoms.position[i] .- sys.atoms.position[k]
#                 rᵢₖ = nearest_mirror!(rᵢₖ, sys.box_sizes_SC)
#                 dist_ik_sq = sum(x -> x^2, rᵢₖ)

#                 if dist_ik_sq < r_cut_sq
#                     nearest_j .= sys.atoms.position[i] .- rᵢⱼ
#                     nearest_k .= sys.atoms.position[i] .- rᵢₖ
        
#                     r_arr .= ustrip.([sys.atoms.position[A]; nearest_j; nearest_k])
#                     block .= H3_exec_ij(r_arr)
#                     IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(j-1) + 1 : D*(j-1) + D] .+= block
#                     IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(i-1) + 1 : D*(i-1) + D] .+= block
#                 end
#             end

#             #When i,j = 1st,3rd atoms in triplet
#             for j in range(1,k-1)
#                 nearest_i = rᵢⱼ .+ sys.atoms.position[j]
#             end

#             #When i,j = 2nd,3rd atoms in triplet
#             for i in range(1,N_atoms)

#             end

#         end
#     end

#     return IFC2

# end