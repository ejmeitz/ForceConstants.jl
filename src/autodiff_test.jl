export second_order_AD_test, third_order_AD_test, fourth_order_AD_test

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
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.box_sizes_SC)

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

    IFC2 = apply_tols!(IFC2,tol)

    return SecondOrderMatrix(IFC2, unit(pot.ϵ / pot.σ^2), tol)

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
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.box_sizes_SC)


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
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.box_sizes_SC)

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

    IFC3 = apply_tols!(IFC3,tol)

    return ThirdOrderMatrix(IFC3, unit(pot.ϵ / pot.σ^3), tol)

end

function fourth_order_AD_test(sys::SuperCellSystem{D}, pot::PairPotential, tol; float_type = Float32) where D
    vars = make_variables(:r, D)
    r_norm = sqrt(sum(x -> x^2, vars))
    pot_symbolic = potential_nounits(pot, r_norm)
    H_symbolic = hessian(pot_symbolic, vars)

    fourth_order_symbolic = reshape(jacobian(vec(jacobian(vec(H_symbolic), vars)), vars),(D,D,D,D))
    FO_exec = make_function(fourth_order_symbolic, vars)

    N_atoms = n_atoms(sys)
    IFC4 = FC_val{float_type,4}[]

    for i in range(1, N_atoms)
        for j in range(i + 1, N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.box_sizes_SC)


            iiij_block = FO_exec(ustrip.(rᵢⱼ))

            for α in range(1,D)
                for β in range(1,D)
                    for γ in range(1,D)
                        for δ in range(1,D)
                            if abs(iiij_block[α,β,γ,δ]) > tol
                                #iiij
                                ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(i-1) + γ; ll = D*(j-1) + δ
                                #& THERE WILL BE DUPLICATES SINCE THIS USES PUSH INSTEAD OF JUST OVERWRITTING DATA!!!
                                set_terms_fourth_order!(χ, ii, jj, kk, ll, float_type(iiij_block[α,β,γ,δ]))
                                #iijj
                                ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(j-1) + γ; ll = D*(j-1) + δ
                                set_terms_fourth_order!(χ, ii, jj, kk, ll, float_type(iiij_block[α,β,γ,δ]))
                            end
                        end
                    end
                end
            end



            rᵢⱼ = sys.atoms.position[j] .- sys.atoms.position[i]
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.box_sizes_SC)

            ijjj_block = FO_exec(ustrip.(rᵢⱼ))

            for α in range(1,D)
                for β in range(1,D)
                    for γ in range(1,D)
                        for δ in range(1,D)
                            if abs(iiij_block[α,β,γ,δ]) > tol
                                #ijjj
                                ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(j-1) + γ; ll = D*(j-1) + δ
                                set_terms_fourth_order!(χ, ii, jj, kk, ll, float_type(ijjj_block[α,β,γ,δ]))
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