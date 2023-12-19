# export fourth_order_AD_test #*export once tests

#*add r_cut 
# function fourth_order(sys::SuperCellSystem{D}, pot::PairPotential,
#     calc::AutoDiffCalculator) where D
#     vars = make_variables(:r, D)
#     r_norm = sqrt(sum(x -> x^2, vars))
#     pot_symbolic = potential_nounits(pot, r_norm)
#     H_symbolic = hessian(pot_symbolic, vars)

#     fourth_order_symbolic = reshape(jacobian(vec(jacobian(vec(H_symbolic), vars)), vars),(D,D,D,D))
#     FO_exec = make_function(fourth_order_symbolic, vars)

#     N_atoms = n_atoms(sys)
#     IFC4 = FC_val{Float64,4}[]

#     for i in range(1, N_atoms)
#         for j in range(i + 1, N_atoms)

#             rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
#             rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)


#             iiij_block = FO_exec(ustrip.(rᵢⱼ))

#             for α in range(1,D)
#                 for β in range(1,D)
#                     for γ in range(1,D)
#                         for δ in range(1,D)
#                             if abs(iiij_block[α,β,γ,δ]) > tol
#                                 #iiij
#                                 ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ; ll = D*(j-1) + δ
#                                 #& THERE WILL BE DUPLICATES SINCE THIS USES PUSH INSTEAD OF JUST OVERWRITTING DATA!!!
#                                 set_terms_fourth_order!(χ, ii, jj, kk, ll, iiij_block[α,β,γ,δ])
#                                 #iijj
#                                 ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(j-1) + γ; ll = D*(j-1) + δ
#                                 set_terms_fourth_order!(χ, ii, jj, kk, ll, iiij_block[α,β,γ,δ])
#                             end
#                         end
#                     end
#                 end
#             end



#             rᵢⱼ = sys.atoms.position[j] .- sys.atoms.position[i]
#             rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)

#             ijjj_block = FO_exec(ustrip.(rᵢⱼ))

#             for α in range(1,D)
#                 for β in range(1,D)
#                     for γ in range(1,D)
#                         for δ in range(1,D)
#                             if abs(iiij_block[α,β,γ,δ]) > tol
#                                 #ijjj
#                                 ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(j-1) + γ; ll = D*(j-1) + δ
#                                 set_terms_fourth_order!(χ, ii, jj, kk, ll, ijjj_block[α,β,γ,δ])
#                             end
#                         end
#                     end
#                 end
#             end

#         end
#     end

#     #Acoustic Sum Rule #&TODO
#     # Threads.@threads for i in range(1, N_atoms) # index of block matrix
#     #     for α in range(1,D)
#     #         for β in range(1,D)
#     #             ii = D*(i-1) + α
#     #             jj = D*(i-1) + β # i == j because we're on diagonal
#     #             IFC2[ii,jj] = -1*sum(IFC2[ii, β:D:end])
#     #         end
#     #     end
#     # end


#     return IFC4

# end
