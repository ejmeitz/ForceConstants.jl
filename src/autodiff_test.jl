export second_order_AD_test, third_order_AD_test, fourth_order_AD_test,energy_loop_mcc

function second_order_AD_test(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    vars = make_variables(:r, D)
    r_norm = sqrt(sum(x -> x^2, vars))
    pot_symbolic = potential_nounits(pot, r_norm)
    H_symbolic = hessian(pot_symbolic, vars)
    H_exec = make_function(H_symbolic, vars)

    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)

    for i in range(1, N_atoms)
        for j in range(i + 1, N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)

            ij_block = H_exec(ustrip.(rᵢⱼ))

            IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(j-1) + 1 : D*(j-1) + D] .= ij_block
            IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(i-1) + 1 : D*(i-1) + D] .= ij_block

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

    IFC2 = apply_tols!(IFC2,tol)

    return DenseForceConstants(IFC2, energy_unit(pot) / length_unit(pot)^2, tol)

end

function second_order_AD_test(sys::SuperCellSystem{D}, pot::StillingerWeberSilicon, tol; r_cut = pot.r_cut) where D

    #2nd order part
    r_ij_vars = make_variables(:rij, D) #rᵢ - rⱼ
    r_ij_norm = sqrt(sum(x -> x^2, r_ij_vars))

    pot2_symbolic = pair_potential_nounits(pot, r_ij_norm)
    H2_symbolic = hessian(pot2_symbolic, r_ij_vars)
    H2_exec = make_function(H2_symbolic, r_ij_vars)


    # Three body part need to do a bit differently
    r_i_vars = make_variables(:ri, D)
    r_j_vars = make_variables(:rj, D)
    r_k_vars = make_variables(:rk, D)

    r_ij_threebody = r_i_vars .- r_j_vars
    r_ik_threebody = r_i_vars .- r_k_vars
    r_ij_norm_threebody = sqrt(sum(x -> x^2, r_ij_threebody))
    r_ik_norm_threebody = sqrt(sum(x -> x^2, r_ik_threebody))

    #3rd order part
    pot3_symbolic = three_body_potential_nounits(pot, r_ij_threebody, r_ik_threebody, r_ij_norm_threebody, r_ik_norm_threebody)
    
    #Can't treat block of IFC as hessian
    H3_symbolic_ij = Matrix{FastDifferentiation.Node}(undef, D, D)
    for a in range(1,D)
        for b in range(a,D)
            H3_symbolic_ij[a,b] = derivative([pot3_symbolic], r_i_vars[a], r_j_vars[b])[1]
            H3_symbolic_ij[b,a] = H3_symbolic_ij[a,b]
        end
    end

    H3_symbolic_jk = Matrix{FastDifferentiation.Node}(undef, D, D)
    for a in range(1,D)
        for b in range(a,D)
            H3_symbolic_jk[a,b] = derivative([pot3_symbolic], r_j_vars[a], r_k_vars[b])[1]
            H3_symbolic_jk[b,a] = H3_symbolic_jk[a,b]
        end
    end

    H3_symbolic_ik = Matrix{FastDifferentiation.Node}(undef, D, D)
    for a in range(1,D)
        for b in range(a,D)
            H3_symbolic_ik[a,b] = derivative([pot3_symbolic], r_i_vars[a], r_k_vars[b])[1]
            H3_symbolic_ik[b,a] = H3_symbolic_ik[a,b]
        end
    end

    H3_exec_ij = make_function(H3_symbolic_ij, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_jk = make_function(H3_symbolic_jk, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_ik = make_function(H3_symbolic_ik, [r_i_vars; r_j_vars; r_k_vars])
    
    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)
    r_cut_sq = r_cut*r_cut
    #Loop over block matricies in IFC2 Matrix
    # Threads.@threads for a in 1:N_atoms
    #     block = zeros(D,D)
    #     for b in 1:N_atoms #& re write to only loop top half of matrix
    #         if a != b #Ignore diagonal, do that with ASR

    #             if a < b
    #                 #2-body part of potential derivative only non-zero on r_ab term
    #                 r_ab = sys.atoms.position[a] .- sys.atoms.position[b]
    #                 r_ab = nearest_mirror!(r_ab, sys.box_sizes_SC)
    #                 dist_ab = norm(r_ab)

    #                 if dist_ab < r_cut
    #                     block .= H2_exec(ustrip.(r_ab))

    #                     IFC2[D*(a-1) + 1 : D*(a-1) + D, D*(b-1) + 1 : D*(b-1) + D] .+= block
    #                     IFC2[D*(b-1) + 1 : D*(b-1) + D, D*(a-1) + 1 : D*(a-1) + D] .+= block
    #                 end
    #             end

    #             #3-body part of potential derivative non-zero in multiple cases

    #             #First term -- a,b = i,j
    #             for k in range(b + 1, N_atoms)
    #                 if k != a #avoid i = k a != b already checked
    #                     r_ab = sys.atoms.position[a] .- sys.atoms.position[b]
    #                     r_ab = nearest_mirror!(r_ab, sys.box_sizes_SC)
    #                     dist_ab = norm(r_ab)

    #                     r_ak = sys.atoms.position[a] .- sys.atoms.position[k]
    #                     r_ak = nearest_mirror!(r_ak, sys.box_sizes_SC)
    #                     dist_ak = norm(r_ak)

    #                     if dist_ab < r_cut && dist_ak < r_cut
    #                         nearest_j = sys.atoms.position[a] - r_ab
    #                         nearest_k = sys.atoms.position[a] - r_ak
    #                         block .= H3_exec_ij(ustrip.([sys.atoms.position[a]; nearest_j; nearest_k]))
    #                         IFC2[D*(a-1) + 1 : D*(a-1) + D, D*(b-1) + 1 : D*(b-1) + D] .+= block
    #                         IFC2[D*(b-1) + 1 : D*(b-1) + D, D*(a-1) + 1 : D*(a-1) + D] .+= block
    #                     end
    #                 end
    #             end

    #             #Second term -- a,b = i,k
    #             # for j in range(1, b-1) #1-b or 1- Natoms
    #             #     if j != a #avoid i = j, a != b already checked
    #             #         r_aj = sys.atoms.position[a] .- sys.atoms.position[j]
    #             #         r_aj = nearest_mirror!(r_aj, sys.box_sizes_SC)
    #             #         dist_aj = norm(r_aj)

    #             #         r_ab = sys.atoms.position[a] .- sys.atoms.position[b]
    #             #         r_ab = nearest_mirror!(r_ab, sys.box_sizes_SC) 
    #             #         dist_ab = norm(r_ab)
                        

    #             #         if dist_aj < r_cut && dist_ab < r_cut
    #             #             nearest_j = sys.atoms.position[a] - r_aj
    #             #             nearest_k = sys.atoms.position[a] - r_ab
    #             #             block .= H3_exec_ik(ustrip.([sys.atoms.position[a]; nearest_j; nearest_k]))
    #             #             IFC2[D*(a-1) + 1 : D*(a-1) + D, D*(b-1) + 1 : D*(b-1) + D] .+= block
    #             #             # IFC2[D*(b-1) + 1 : D*(b-1) + D, D*(a-1) + 1 : D*(a-1) + D] .+= block
    #             #         end
    #             #     end
    #             # end

    #             #Third term -- a,b = j,k
    #             if a < b
    #                 for i in range(1, N_atoms)
    #                     if i != b && i != a
    #                         r_ia = sys.atoms.position[i] .- sys.atoms.position[a]
    #                         r_ia = nearest_mirror!(r_ia, sys.box_sizes_SC)
    #                         dist_ia = norm(r_ia)
    
    #                         r_ib = sys.atoms.position[i] .- sys.atoms.position[b]
    #                         r_ib = nearest_mirror!(r_ib, sys.box_sizes_SC)
    #                         dist_ib = norm(r_ib)
                            
    
    #                         if dist_ia < r_cut && dist_ib < r_cut
    #                             nearest_j = sys.atoms.position[i] - r_ia
    #                             nearest_k = sys.atoms.position[i] - r_ib
    #                             block .= H3_exec_jk(ustrip.([sys.atoms.position[i]; nearest_j; nearest_k]))
    #                             IFC2[D*(a-1) + 1 : D*(a-1) + D, D*(b-1) + 1 : D*(b-1) + D] .+= block
    #                             # IFC2[D*(b-1) + 1 : D*(b-1) + D, D*(a-1) + 1 : D*(a-1) + D] .+= block
    #                         end
    #                     end
    #                 end
    #             end

    #         end
    #     end
    # end

    # function print_derivs(N_atoms)
    #     for a in 1:N_atoms
    #         for b in 1:N_atoms
    #             if a != b
    #                 if a < b
    #                     println("Partial Φ₂(r$(a)$(b))")
    #                 end

    #                 for k in range(b + 1, N_atoms)
    #                     if k != a #avoid i = k a != b already checked
    #                         println("Partial Φ₃(r$(a)$(b), r$(a)$(k)) from ij")
    #                     end
    #                 end

    #                 # for j in range(1, b-1) #1-b or 1- Natoms
    #                 #     if j != a #avoid i = j, a != b already checked
    #                 #         println("Partial Φ₃(r$(a)$(j), r$(a)$(b)) from ik")
    #                 #     end
    #                 # end

    #                 # if a < b
    #                 #     for i in range(1, N_atoms)
    #                 #         if i != b && i != a
    #                 #             println("Partial Φ₃(r$(i)$(a), r$(i)$(b)) from jk")
    #                 #         end
    #                 #     end
    #                 # end
    #             end
    #         end
    #     end
    # end

   

    #Loop Atomic Interactions and Add their contribution to various derivatives
    Threads.@threads for A in range(1,N_atoms)
        block = zeros(D,D)
        r_ab = similar(sys.atoms.position[1])
        rᵢₖ = similar(sys.atoms.position[1])
        nearest_j = similar(sys.atoms.position[1])
        nearest_k = similar(sys.atoms.position[1])
        for B in range(1, N_atoms)
            if A != B
                #Two body term
                r_ab .= sys.atoms.position[A] .- sys.atoms.position[B]
                nearest_mirror!(r_ab, sys.box_sizes_SC)
                dist_ab_sq = sum(x -> x^2, r_ab)

                if dist_ab_sq < r_cut_sq

                    if A < B
                        block .= H2_exec(ustrip.(r_ab))

                        #*Minus sign needed???
                        IFC2[D*(A-1) + 1 : D*(A-1) + D, D*(B-1) + 1 : D*(B-1) + D] .-= block
                        IFC2[D*(B-1) + 1 : D*(B-1) + D, D*(A-1) + 1 : D*(A-1) + D] .-= block
                    end

                    #Three body terms:
                    for k in range(B+1, N_atoms)
                        if A != k
                            rᵢₖ .= sys.atoms.position[A] .- sys.atoms.position[k]
                            nearest_mirror!(rᵢₖ, sys.box_sizes_SC)
                            dist_ik_sq = sum(x -> x^2, rᵢₖ)
                            # dist_ik_sq = sum(x -> x^2, rᵢₖ)
                            
                            if dist_ik_sq < r_cut_sq
                                #ϕ₃(rᵢⱼ, rᵢₖ) contributes to ∂ri∂rj, ∂ri∂rk, ∂rj∂rk (and the symmetric versions)
                                nearest_j .= sys.atoms.position[A] - r_ab
                                nearest_k .= sys.atoms.position[A] - rᵢₖ
                                #ij terms
                                block .= H3_exec_ij(ustrip.([sys.atoms.position[A]; nearest_j; nearest_k]))
                                IFC2[D*(A-1) + 1 : D*(A-1) + D, D*(B-1) + 1 : D*(B-1) + D] .+= block
                                
                                #*Why does removing tehse do nothing??
                                #ik terms
                                block .= H3_exec_ik(ustrip.([sys.atoms.position[A]; nearest_j; nearest_k]))
                                IFC2[D*(A-1) + 1 : D*(A-1) + D, D*(k-1) + 1 : D*(k-1) + D] .+= block
                                 
                                #jk terms
                                block .= H3_exec_jk(ustrip.([sys.atoms.position[A]; nearest_j; nearest_k]))
                                IFC2[D*(B-1) + 1 : D*(B-1) + D, D*(k-1) + 1 : D*(k-1) + D] .+= block

                                #*do these need symmetries?
                                #* right now only populates top half of matrix cause k > j

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

    IFC2 = apply_tols!(IFC2,tol)

    return DenseForceConstants(IFC2, energy_unit(pot) / length_unit(pot)^2, tol)

end

function third_order_AD_test(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    vars = make_variables(:r, D)
    r_norm = sqrt(sum(x -> x^2, vars))
    pot_symbolic = potential_nounits(pot, r_norm)
    H_symbolic = hessian(pot_symbolic, vars)

    third_order_symbolic = reshape(jacobian(vec(H_symbolic), vars),(D,D,D))
    TO_exec = make_function(third_order_symbolic, vars)

    N_atoms = n_atoms(sys)
    IFC3 = zeros(D*N_atoms,D*N_atoms,D*N_atoms)

    for i in range(1, N_atoms)
        for j in range(i + 1, N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)


            iij_block = TO_exec(ustrip.(rᵢⱼ))

            #i,i,j terms
            IFC3[D*(i-1) + 1 : D*(i-1) + D,
                 D*(i-1) + 1 : D*(i-1) + D,
                 D*(j-1) + 1 : D*(j-1) + D] .= iij_block

            IFC3[D*(i-1) + 1 : D*(i-1) + D,
                 D*(j-1) + 1 : D*(j-1) + D,
                 D*(i-1) + 1 : D*(i-1) + D] .= iij_block
            
            IFC3[D*(j-1) + 1 : D*(j-1) + D,
                 D*(i-1) + 1 : D*(i-1) + D,
                 D*(i-1) + 1 : D*(i-1) + D] .= iij_block


            rᵢⱼ = sys.atoms.position[j] .- sys.atoms.position[i]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)

            ijj_block = TO_exec(ustrip.(rᵢⱼ))

            #ijj terms
            IFC3[D*(i-1) + 1 : D*(i-1) + D,
                 D*(j-1) + 1 : D*(j-1) + D,
                 D*(j-1) + 1 : D*(j-1) + D] .= ijj_block

            IFC3[D*(j-1) + 1 : D*(j-1) + D,
                 D*(j-1) + 1 : D*(j-1) + D,
                 D*(i-1) + 1 : D*(i-1) + D] .= ijj_block
            
            IFC3[D*(j-1) + 1 : D*(j-1) + D,
                 D*(i-1) + 1 : D*(i-1) + D,
                 D*(j-1) + 1 : D*(j-1) + D] .= ijj_block

        end
    end

    #Acoustic Sum Rule
    ASR!(IFC3, N_atoms, D)

    IFC3 = apply_tols!(IFC3,tol)

    return DenseForceConstants(IFC3, energy_unit(pot) / length_unit(pot)^3, tol)

end

function fourth_order_AD_test(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    vars = make_variables(:r, D)
    r_norm = sqrt(sum(x -> x^2, vars))
    pot_symbolic = potential_nounits(pot, r_norm)
    H_symbolic = hessian(pot_symbolic, vars)

    fourth_order_symbolic = reshape(jacobian(vec(jacobian(vec(H_symbolic), vars)), vars),(D,D,D,D))
    FO_exec = make_function(fourth_order_symbolic, vars)

    N_atoms = n_atoms(sys)
    IFC4 = FC_val{Float64,4}[]

    for i in range(1, N_atoms)
        for j in range(i + 1, N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)


            iiij_block = FO_exec(ustrip.(rᵢⱼ))

            for α in range(1,D)
                for β in range(1,D)
                    for γ in range(1,D)
                        for δ in range(1,D)
                            if abs(iiij_block[α,β,γ,δ]) > tol
                                #iiij
                                ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ; ll = D*(j-1) + δ
                                #& THERE WILL BE DUPLICATES SINCE THIS USES PUSH INSTEAD OF JUST OVERWRITTING DATA!!!
                                set_terms_fourth_order!(χ, ii, jj, kk, ll, iiij_block[α,β,γ,δ])
                                #iijj
                                ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(j-1) + γ; ll = D*(j-1) + δ
                                set_terms_fourth_order!(χ, ii, jj, kk, ll, iiij_block[α,β,γ,δ])
                            end
                        end
                    end
                end
            end



            rᵢⱼ = sys.atoms.position[j] .- sys.atoms.position[i]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)

            ijjj_block = FO_exec(ustrip.(rᵢⱼ))

            for α in range(1,D)
                for β in range(1,D)
                    for γ in range(1,D)
                        for δ in range(1,D)
                            if abs(iiij_block[α,β,γ,δ]) > tol
                                #ijjj
                                ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(j-1) + γ; ll = D*(j-1) + δ
                                set_terms_fourth_order!(χ, ii, jj, kk, ll, ijjj_block[α,β,γ,δ])
                            end
                        end
                    end
                end
            end

        end
    end

    #Acoustic Sum Rule #&TODO
    # Threads.@threads for i in range(1, N_atoms) # index of block matrix
    #     for α in range(1,D)
    #         for β in range(1,D)
    #             ii = D*(i-1) + α
    #             jj = D*(i-1) + β # i == j because we're on diagonal
    #             IFC2[ii,jj] = -1*sum(IFC2[ii, β:D:end])
    #         end
    #     end
    # end


    return IFC4

end


function energy_loop_mcc(sys::SuperCellSystem{D}, pot::PairPotential, q, phi, mass_no_units) where D

    inv_sqrt_masses = 1.0./(sqrt.(mass_no_units))
    N_atoms = n_atoms(sys)
    N_modes = D*N_atoms

    posns = [inv_sqrt_masses[(i%3) + 1]*dot(q, phi[i+1,:]) for i in 0:N_modes-1]

    U_total = 0.0
    r_cut = ustrip(pot.r_cut)
    box_sizes = ustrip.(sys.box_sizes_SC)

    for i in range(1,N_atoms)
        for j in range(i+1, N_atoms)             
            r = posns[i:i+D-1] .- posns[j:j+D-1]
            nearest_mirror!(r, box_sizes)
            dist = norm(r)
            if dist < r_cut
                U_total += potential_nounits(pot,dist)
            end
        end
    end

    return U_total

end

#q_vars = make_variable(:q, N_modes)
#r_vars = zeros(eltype(q_vars), length(q_vars))
#for i in 1:N_modes
#   r_vars[i] = dot(q_vars, phi[i,:])
#end

# divide r_vars by mass_sqrt

# pass r_vars into energy_loop2

# function energy_loop2(pot::PairPotential, posns, r_cut, box_sizes, N_atoms)

#     U_total = 0.0
#     box_sizes = ustrip.(box_sizes)

#     for i in range(1,N_atoms)
#         for j in range(i+1, N_atoms)             
#             r = posns[3*(i-1) + 1: 3*i] .- posns[3*(j-1) + 1: 3*j ]
#             nearest_mirror_AD!(r, box_sizes)
#             dist = sqrt(sum(x-> x*x, r))
#             # if dist < r_cut
#             U_total += potential_nounits(pot, dist)
#             # end
#         end
#     end

#     return U_total

# end

# #& not quite right but I just wanna test the AD
# function nearest_mirror_AD!(r_ij, box_sizes)
#     r_ij .+= sign.(r_ij .- box_sizes./2) .* box_sizes
#     return r_ij
# end