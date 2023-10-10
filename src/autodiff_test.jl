export second_order_AD_test, third_order_AD_test

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

    third_order_symbolic = zeros(eltype(H_symbolic), (D,D,D))
    #Build third order tensor from hessian
    for i in range(1,D)
        for j in range(1,D)
            for k in range(1,D)
                third_order_symbolic[i,j,k] = derivative([H_symbolic[i,j];],vars[k])
            end
        end
    end
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
                 D*(j-1) + 1 : D*(j-1) + D] .= ijk_block


            rᵢⱼ = sys.atoms.position[j] .- sys.atoms.position[i]
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.box_sizes_SC)

            ijj_block = TO_exec(ustrip.(rᵢⱼ))

            #ijj terms
            IFC3[D*(i-1) + 1 : D*(i-1) + D,
                 D*(j-1) + 1 : D*(j-1) + D,
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