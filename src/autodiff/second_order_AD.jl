function second_order!(IFC2::Matrix{T}, sys::SuperCellSystem{D}, pot::PairPotential,
      calc::AutoDiffCalculator) where {T,D}

    @assert calc.r_cut <= pot.r_cut "Calculator r_cut must be less than potential r_cut"  

    vars = make_variables(:r, D)
    r_norm = sqrt(sum(x -> x^2, vars))
    pot_symbolic = potential_nounits(pot, r_norm)
    H_symbolic = hessian(pot_symbolic, vars)
    H_exec = make_function(H_symbolic, vars)
    r_cut_sq = calc.r_cut*calc.r_cut

    N_atoms = n_atoms(sys)

    for i in range(1, N_atoms)
        for j in range(i + 1, N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
            dist_ij_sq = sum(x -> x^2, rᵢⱼ)

            if dist_ij_sq < r_cut_sq
                ij_block = -H_exec(ustrip.(rᵢⱼ))

                IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(j-1) + 1 : D*(j-1) + D] .= ij_block
                IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(i-1) + 1 : D*(i-1) + D] .= ij_block
            end
        end
    end

    ASR!(IFC2, N_atoms, D)

    return IFC2

end

function second_order!(IFC2::Matrix{T}, sys::SuperCellSystem{D}, pot::StillingerWeberSilicon,
     calc::AutoDiffCalculator) where {D,T}

    @assert calc.r_cut <= pot.r_cut "For SW silicon force constant 
        cutoff must be less than potential cutoff"

    H2_exec, H3_exec_ij, H3_exec_ik, H3_exec_jk = 
        three_body_second_derivs(pot, D)
   
    N_atoms = n_atoms(sys)
    r_cut_sq = calc.r_cut*calc.r_cut  

    # ij_locks = Array{ReentrantLock}(undef, D*N_atoms, D*N_atoms)
    # fill!(ij_locks, ReentrantLock())

    # jk_locks = Array{ReentrantLock}(undef, D*N_atoms, D*N_atoms)
    # fill!(jk_locks, ReentrantLock())

    #Loop Atomic Interactions and Add their contribution to various derivatives
    for i in range(1,N_atoms)
        block = zeros(D,D)
        block2 = zeros(D,D)
        rᵢⱼ = similar(sys.atoms.position[1])
        rᵢₖ = similar(sys.atoms.position[1])
        nearest_j = similar(sys.atoms.position[1])
        nearest_k = similar(sys.atoms.position[1])
        r_arr = Vector{Float64}(undef, D*D)
        for j in range(1, N_atoms)
            if i != j
                rᵢⱼ .= sys.atoms.position[i] .- sys.atoms.position[j]
                nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
                dist_ab_sq = sum(x -> x^2, rᵢⱼ)

                if dist_ab_sq < r_cut_sq
                    #Two body term only contributes to ij and ji blocks
                    if j > i
                        block .= -H2_exec(ustrip.(rᵢⱼ))
                        IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(j-1) + 1 : D*(j-1) + D] .+= block
                        IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(i-1) + 1 : D*(i-1) + D] .+= block
                    end

                    #Three body terms between atoms i,j,k
                    for k in range(j+1, N_atoms)
                        if i != k
                            rᵢₖ .= sys.atoms.position[i] .- sys.atoms.position[k]
                            nearest_mirror!(rᵢₖ, sys.box_sizes_SC)
                            dist_ik_sq = sum(x -> x^2, rᵢₖ)
                            
                            if dist_ik_sq < r_cut_sq
                                nearest_j .= sys.atoms.position[i] .- rᵢⱼ
                                nearest_k .= sys.atoms.position[i] .- rᵢₖ
                 
                                #contribution to ij derivative block
                                r_arr .= ustrip.([sys.atoms.position[i]; nearest_j; nearest_k])
                                block .= H3_exec_ij(r_arr)
                                # lock(ij_locks[D*(i-1) + 1, D*(j-1) + 1])
                                IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(j-1) + 1 : D*(j-1) + D] .+= block
                                IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(i-1) + 1 : D*(i-1) + D] .+= permutedims!(block2, block, (2,1))
                                # unlock(ij_locks[D*(i-1) + 1, D*(j-1) + 1])

                                #contribution to ik derivative block
                                block .= H3_exec_ik(r_arr)
                                IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(k-1) + 1 : D*(k-1) + D] .+= block
                                IFC2[D*(k-1) + 1 : D*(k-1) + D, D*(i-1) + 1 : D*(i-1) + D] .+= permutedims!(block2, block, (2,1))

                                block .= H3_exec_jk(r_arr)
                                # lock(jk_locks[D*(j-1) + 1, D*(k-1) + 1])
                                IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(k-1) + 1 : D*(k-1) + D] .+= block
                                IFC2[D*(k-1) + 1 : D*(k-1) + D, D*(j-1) + 1 : D*(j-1) + D] .+= permutedims!(block2, block, (2,1))
                                # unlock(jk_locks[D*(j-1) + 1, D*(k-1) + 1])

                            end
                        end
                    end
                end
            end

        end
    end

    #Acoustic Sum Rule
    ASR!(IFC2, N_atoms, D)
    
    return IFC2

end

