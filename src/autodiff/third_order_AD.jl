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

    for i in range(1, N_atoms)
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


               #& THIS EQUIVALENT TO ABOVE??
               rᵢⱼ = sys.atoms.position[j] .- sys.atoms.position[i]
               rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)

               ijj_block = -TO_exec(ustrip.(rᵢⱼ))

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

     H2_exec, H3_exic_iij, H3_exic_iik,
     H3_exic_ijj, H3_exic_ijk, H3_exic_ikk = 
        three_body_third_derivs(pot, D)
   
    N_atoms = n_atoms(sys)
    IFC3 = zeros(D*N_atoms, D*N_atoms, D*N_atoms)
    r_cut_sq = calc.r_cut*calc.r_cut  

    #Loop Atomic Interactions and Add their contribution to various derivatives
    Threads.@threads for A in range(1,N_atoms)
        block = zeros(D,D,D)
        r_ab = similar(sys.atoms.position[1])
        r_ba = similar(sys.atoms.position[1])
        rᵢₖ = similar(sys.atoms.position[1])
        nearest_j = similar(sys.atoms.position[1])
        nearest_k = similar(sys.atoms.position[1])
        r_arr = Vector{Float64}(undef, D*D)

        i_rng = D*(A-1) + 1 : D*(A-1) + D
        for B in range(1, N_atoms)
            if A != B
                j_rng = D*(B-1) + 1 : D*(B-1) + D #*this allocates new range every time
                #Two body term
                r_ab .= sys.atoms.position[A] .- sys.atoms.position[B]
                nearest_mirror!(r_ab, sys.box_sizes_SC)
                dist_ab_sq = sum(x -> x^2, r_ab)

                if dist_ab_sq < r_cut_sq

                    if A < B
                        block .= -H2_exec(ustrip.(r_ab))
                        #iij
                        set_third_order_terms!(IFC3, i_rng, j_rng, block)
                        
                        #ijj
                        r_ba .= sys.atoms.position[B] .- sys.atoms.position[A]
                        nearest_mirror!(r_ba, sys.box_sizes_SC)

                        block .= -H2_exec(ustrip.(r_ba))
                        set_third_order_terms!(IFC3, j_rng, i_rng, block)
                    end

                    #Three body terms:
                    for k in range(B+1, N_atoms)
                        if A != k
                            k_rng = D*(k-1) + 1 : D*(k-1) + D
                            rᵢₖ .= sys.atoms.position[A] .- sys.atoms.position[k]
                            nearest_mirror!(rᵢₖ, sys.box_sizes_SC)
                            dist_ik_sq = sum(x -> x^2, rᵢₖ)
                            
                            if dist_ik_sq < r_cut_sq
                                nearest_j .= sys.atoms.position[A] .- r_ab
                                nearest_k .= sys.atoms.position[A] .- rᵢₖ
                 
                                #contribution to ij derivative block
                                r_arr .= ustrip.([sys.atoms.position[A]; nearest_j; nearest_k])

                                
                                block .= H3_exic_iij(r_arr)
                                set_third_order_terms!(IFC3, i_rng, j_rng, block)
                                
                                #&TODO figure out how to to the jji term elegantly
                                # #ijj
                                # block .= -H2_exec(ustrip.(r_ba))
                                # set_third_order_terms!(IFC3, j_rng, i_rng, block)


                                
                                #contribution to ik derivative block
                                block .= H3_exec_ik(r_arr)
                                IFC2[D*(A-1) + 1 : D*(A-1) + D, D*(k-1) + 1 : D*(k-1) + D] .+= block
                                IFC2[D*(k-1) + 1 : D*(k-1) + D, D*(A-1) + 1 : D*(A-1) + D] .+= block
                                 
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