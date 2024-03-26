function third_order!(IFC3::Array{T,3}, sys::SuperCellSystem{D}, pot::PairPotential,
      calc::AutoDiffCalculator) where {T,D}

    vars = make_variables(:r, D)
    r_norm = sqrt(sum(x -> x^2, vars))
    pot_symbolic = potential_nounits(pot, r_norm)
    H_symbolic = hessian(pot_symbolic, vars)

    third_order_symbolic = reshape(jacobian(vec(H_symbolic), vars),(D,D,D))
    TO_exec = make_function(third_order_symbolic, vars)

    N_atoms = n_atoms(sys)
    r_cut_sq = calc.r_cut*calc.r_cut

    Threads.@threads for i in range(1, N_atoms)
        for j in range(i + 1, N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
            dist_ij_sq = sum(x -> x^2, rᵢⱼ)

            if dist_ij_sq < r_cut_sq
                
               iij_block = -TO_exec(ustrip.(rᵢⱼ))

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


               #2 js --> negate result again
               ijj_block = -iij_block #deriv wrt to r_ij so can re-use above

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
    end

    #Acoustic Sum Rule
    ASR!(IFC3, N_atoms, D)

    return IFC3

end


function third_order!(IFC3::Array{T,3}, sys::SuperCellSystem{D}, pot::StillingerWeberSilicon,
     calc::AutoDiffCalculator) where {T,D}

     @assert calc.r_cut <= pot.r_cut "For AutoDiff SW silicon force constant 
        cutoff must be less than potential cutoff"

    N_atoms = n_atoms(sys)
    r_cut_sq = calc.r_cut*calc.r_cut  

    H3_exec_iik = calc_H3_exec_iik(pot, D) #FastDiff breaks on this one so use Symbolics.jl (slow)

    block = zeros(D,D,D)
    block2 = zeros(D,D,D)
    rᵢⱼ = similar(sys.atoms.position[1])
    rᵢₖ = similar(sys.atoms.position[1])
    nearest_j = similar(sys.atoms.position[1])
    nearest_k = similar(sys.atoms.position[1])
    r_arr = Vector{Float64}(undef, D*D)

    #Loop Atomic Interactions and accumulate their contribution to various derivatives
    for i in range(1,N_atoms) #&is this safe to parallelize?
        i_rng = D*(i-1) + 1 : D*(i-1) + D
        for j in range(1, N_atoms)
            if i != j
                j_rng = D*(j-1) + 1 : D*(j-1) + D
                #Two body term
                rᵢⱼ .= sys.atoms.position[i] .- sys.atoms.position[j]
                nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
                dist_ij_sq = sum(x -> x^2, rᵢⱼ)

                if dist_ij_sq < r_cut_sq

                    if j > i #two body derivatives wrt r_ija, r_ijb, r_ijk
                        block .= H3_exec_ij(ustrip.(rᵢⱼ)) 

                        #iij
                        set_third_order_terms!(IFC3, i_rng, j_rng, -block)
                        
                        #ijj
                        set_third_order_terms!(IFC3, j_rng, i_rng, block)
                    end

                    #Three body terms:
                    for k in range(j+1, N_atoms)
                        if i != k
                            k_rng = D*(k-1) + 1 : D*(k-1) + D
                            rᵢₖ .= sys.atoms.position[i] .- sys.atoms.position[k]
                            nearest_mirror!(rᵢₖ, sys.box_sizes_SC)
                            dist_ik_sq = sum(x -> x^2, rᵢₖ)
                            
                            if dist_ik_sq < r_cut_sq #three body derivatives wrt r_ia, r_ib, r_jk
                                nearest_j .= sys.atoms.position[i] .- rᵢⱼ
                                nearest_k .= sys.atoms.position[i] .- rᵢₖ

                                update_r_arr!(r_arr, sys.atoms.position[i], nearest_j, nearest_k, D)
                                # r_arr .= ustrip.([sys.atoms.position[i]; nearest_j; nearest_k])

                                # #ij #*These are not thread safe since i and j can be same perm (e.g. (1,2), (2,1))
                                block .= H3_exec_iij(r_arr)
                                # lock(ij_locks[D*(i-1) + 1, D*(j-1) + 1])
                                set_third_order_terms!(IFC3, i_rng, j_rng, block, block2)

                                block .= H3_exec_ijj(r_arr)
                                set_third_order_terms_alt!(IFC3, j_rng, i_rng, block, block2)
                                # unlock(ij_locks[D*(i-1) + 1, D*(j-1) + 1])

                                # #ik
                                # block .= H3_exec_iik(r_arr) #* CAUSES ASYMMETRY
                                # set_third_order_terms!(IFC3, i_rng, k_rng, block, block2)

                                H3_exec_symb!(block, H3_exec_iik, r_arr)
                                set_third_order_terms!(IFC3, i_rng, k_rng, block, block2)

                                block .= H3_exec_ikk(r_arr)
                                set_third_order_terms_alt!(IFC3, k_rng, i_rng, block, block2)

                                # #jk #*These are not thread safe since there is no i
                                block .= H3_exec_jjk(r_arr)
                                set_third_order_terms!(IFC3, j_rng, k_rng, block, block2)

                                block .= H3_exec_jkk(r_arr)
                                set_third_order_terms_alt!(IFC3, k_rng, j_rng, block, block2)

                                #ijk
                                block .= H3_exec_ijk(r_arr)
                                set_third_order_terms!(IFC3, i_rng, j_rng, k_rng, block, block2)
                            end
                        end
                    end
                end
            end

        end
    end

    IFC3 = ASR!(IFC3, N_atoms, D)

    return IFC3

end

function update_r_arr!(r_arr, r_i, r_j, r_k, D)
    r_arr[1:D] .= ustrip(r_i)
    r_arr[D+1:2*D] .= ustrip(r_j)
    r_arr[2*D+1:3*D] .= ustrip(r_k)
    return r_arr
end


function set_third_order_terms!(arr, rng1::UnitRange, rng2::UnitRange, block)

    arr[rng1,rng1,rng2] .+= block
    arr[rng1,rng2,rng1] .+= block
    arr[rng2,rng1,rng1] .+= block

    return arr

end

#aab terms
function set_third_order_terms!(arr, rng1::UnitRange, rng2::UnitRange, block, block2)

    arr[rng1,rng1,rng2] .+= block
    arr[rng1,rng2,rng1] .+= permutedims!(block2, block, (1,3,2))
    arr[rng2,rng1,rng1] .+= permutedims!(block2, block, (3,2,1))

    return arr

end

#abb terms
function set_third_order_terms_alt!(arr, rng1::UnitRange, rng2::UnitRange, block, block2)

    arr[rng2,rng1,rng1] .+= block
    arr[rng1,rng2,rng1] .+= permutedims!(block2, block, (2,1,3))
    arr[rng1,rng1,rng2] .+= permutedims!(block2, block, (3,2,1))

    return arr

end

function set_third_order_terms!(arr, rng1::UnitRange, rng2::UnitRange,
    rng3::UnitRange, block, block2)

   arr[rng1,rng2,rng3] .+= block
   arr[rng1,rng3,rng2] .+= permutedims!(block2, block, (1,3,2)) 
   arr[rng2,rng1,rng3] .+= permutedims!(block2, block, (2,1,3))
   arr[rng2,rng3,rng1] .+= permutedims!(block2, block, (2,3,1))
   arr[rng3,rng1,rng2] .+= permutedims!(block2, block, (3,1,2))
   arr[rng3,rng2,rng1] .+= permutedims!(block2, block, (3,2,1))

   return arr

end

# #Useful for debugging the rotations in 3D
# function make_pattern(a,b,c)
#     out = Array{Vector{String}}(undef, (3,3,3))
#     for (i,x) in enumerate([:x,:y,:z])
#         for (j,y) in enumerate([:x,:y,:z])
#             for (k,z) in enumerate([:x,:y,:z])
#                 out[i,j,k] = sort(["r_$(a)$(x)", "r_$(b)$(y)", "r_$(c)$(z)"])
#             end
#         end
#     end
#     return out
# end

# # Useful for debugging symmetry issues
# function check_sym(arr, k, N_at, D)
#     for i in 1:N_at
#         for j in i+1:N_at
#             for α in 1:D
#                 for β in 1:D
#                      if !isapprox(arr[D*(i-1) + α, D*(j-1) + β, k], arr[D*(j-1) + β, D*(i-1) + α, k])
#                         println("i: $i, j: $j, k: $(((k-1) ÷ 3) + 1) ")
#                     end
#                 end
#             end
#         end
#     end
# end

# function check_wrong(arr, arr2, k, N_at, D)
#     for i in 1:N_at
#         for j in i+1:N_at
#             for α in 1:D
#                 for β in 1:D
#                     v1 = arr[D*(i-1) + α, D*(j-1) + β, k]
#                     v2 = arr2[D*(j-1) + β, D*(i-1) + α, k]
#                     if !isapprox(v1, v2, atol = 1e-4)
#                         println("i: $i, j: $j, k: $(((k-1) ÷ 3) + 1), err:  $(abs(v2 - v1))")
#                     end
#                 end
#             end
#         end
#     end
# end

# function count_wrong(arr, arr2, tol = 1e-4)
#     return sum(isapprox.(arr, arr2, atol = tol) .== false)
# end

# function wrong_dist_hist(arr, arr2, sys)

#     N_atoms = n_atoms(sys)
#     D = 3

#     for i in 1:N_atoms
#         for j in 1:N_atoms
#             for k in [1]
#                 for α in 1:D
#                     for β in 1:D
#                         for γ in 1:D
#                             v1 = arr[D*(i-1) + α, D*(j-1) + β, D*(k-1) + γ]
#                             v2 = arr2[D*(j-1) + β, D*(i-1) + α, D*(k-1) + γ]
#                             if !isapprox(v1, v2, atol = 1e-4)
#                                 # println("i: $i, j: $j, k: $(((k-1) ÷ 3) + 1), err:  $(abs(v2 - v1))")
                                
#                                 r_ij = sys.atoms.position[i] .- sys.atoms.position[j]
#                                 r_ik = sys.atoms.position[i] .- sys.atoms.position[k]
#                                 dist_12 = round(ustrip(norm(nearest_mirror!(r_ij, sys.box_sizes_SC))), digits = 4)
#                                 dist_13 = round(ustrip(norm(nearest_mirror!(r_ik, sys.box_sizes_SC))), digits = 4)
#                                 println("r_$i,$j: $(dist_12), r_$i,$k: $(dist_13)")
#                             end
#                         end
#                     end
#                 end
#             end
#         end
#     end

# end