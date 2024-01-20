export third_order

function third_order(sys::SuperCellSystem{D}, pot::PairPotential,
      calc::AutoDiffCalculator) where D

    vars = make_variables(:r, D)
    r_norm = sqrt(sum(x -> x^2, vars))
    pot_symbolic = potential_nounits(pot, r_norm)
    H_symbolic = hessian(pot_symbolic, vars)

    third_order_symbolic = reshape(jacobian(vec(H_symbolic), vars),(D,D,D))
    TO_exec = make_function(third_order_symbolic, vars)

    N_atoms = n_atoms(sys)
    IFC3 = zeros(D*N_atoms,D*N_atoms,D*N_atoms)
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


               #2 js --> negatives cancel out
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
    end

    #Acoustic Sum Rule
    ASR!(IFC3, N_atoms, D)

    IFC3 = apply_tols!(IFC3, calc.tol)

    return DenseForceConstants(IFC3, energy_unit(pot) / length_unit(pot)^3, calc.tol)

end


function third_order(sys::SuperCellSystem{D}, pot::StillingerWeberSilicon,
     calc::AutoDiffCalculator) where D

     @assert calc.r_cut <= pot.r_cut "For SW silicon force constant 
        cutoff must be less than potential cutoff"

     H3_exec_ij, H3_exec_iij, H3_exec_iik, H3_exec_ijj, H3_exec_ijk, H3_exec_ikk = 
        three_body_third_derivs(pot, D)
   
    N_atoms = n_atoms(sys)
    IFC3 = zeros(D*N_atoms, D*N_atoms, D*N_atoms)
    r_cut_sq = calc.r_cut*calc.r_cut  

    #Loop Atomic Interactions and Add their contribution to various derivatives
    Threads.@threads for i in range(1,N_atoms) #&is this safe to parallelize?
        block = zeros(D,D,D)
        rᵢⱼ = similar(sys.atoms.position[1])
        rᵢₖ = similar(sys.atoms.position[1])
        nearest_j = similar(sys.atoms.position[1])
        nearest_k = similar(sys.atoms.position[1])
        r_arr = Vector{Float64}(undef, D*D)

        i_rng = D*(i-1) + 1 : D*(i-1) + D
        for j in range(1, N_atoms)
            if i != j
                j_rng = D*(j-1) + 1 : D*(j-1) + D #*this allocates new range every time
                #Two body term
                rᵢⱼ .= sys.atoms.position[i] .- sys.atoms.position[j]
                nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
                dist_ij_sq = sum(x -> x^2, rᵢⱼ)

                if dist_ij_sq < r_cut_sq

                    if j > i
                        #iij
                        block .= -H3_exec_ij(ustrip.(rᵢⱼ))
                        set_third_order_terms!(IFC3, i_rng, j_rng, block)
                        
                        #ijj
                        block .= H3_exec_ij(ustrip.(rᵢⱼ))
                        set_third_order_terms!(IFC3, j_rng, i_rng, block)
                    end

                    # #Three body terms:
                    for k in range(j+1, N_atoms)
                        if i != k
                            k_rng = D*(k-1) + 1 : D*(k-1) + D
                            rᵢₖ .= sys.atoms.position[i] .- sys.atoms.position[k]
                            nearest_mirror!(rᵢₖ, sys.box_sizes_SC)
                            dist_ik_sq = sum(x -> x^2, rᵢₖ)
                            
                            if dist_ik_sq < r_cut_sq
                                nearest_j .= sys.atoms.position[i] .- rᵢⱼ
                                nearest_k .= sys.atoms.position[i] .- rᵢₖ
                                r_arr .= ustrip.([sys.atoms.position[i]; nearest_j; nearest_k])

                                block .= H3_exec_iij(r_arr)
                                set_third_order_terms!(IFC3, i_rng, j_rng, block)

                                block .= H3_exec_iik(r_arr)
                                set_third_order_terms!(IFC3, i_rng, k_rng, block)

                                block .= H3_exec_ijk(r_arr)
                                set_third_order_terms!(IFC3, i_rng, j_rng, k_rng, block)
                                
                                block .= H3_exec_ijj(r_arr)
                                set_third_order_terms!(IFC3, j_rng, i_rng, block)

                                block .= H3_exec_ikk(r_arr)
                                set_third_order_terms!(IFC3, k_rng, i_rng, block)
                            end
                        end
                    end
                end
            end

        end
    end

    IFC3 = ASR!(IFC3, N_atoms, D)

    IFC3 = apply_tols!(IFC3, calc.tol)

    return DenseForceConstants(IFC3, energy_unit(pot) / length_unit(pot)^3, calc.tol)

end


function set_third_order_terms!(arr, rng1::UnitRange, rng2::UnitRange, block)

    arr[rng1,rng1,rng2] .+= block
    arr[rng1,rng2,rng1] .+= block
    arr[rng2,rng1,rng1] .+= block

    return arr

end

function set_third_order_terms!(arr, rng1::UnitRange, rng2::UnitRange,
     rng3::UnitRange, block)

    arr[rng1,rng1,rng3] .+= block
    arr[rng1,rng3,rng2] .+= block
    arr[rng2,rng1,rng3] .+= block
    arr[rng2,rng3,rng1] .+= block
    arr[rng3,rng1,rng2] .+= block
    arr[rng3,rng2,rng1] .+= block

    return arr

end