export third_order

function third_order(sys::SuperCellSystem{D}, pot::PairPotential,
      tol, calc::AutoDiffCalculator) where D

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
    end

    #Acoustic Sum Rule
    ASR!(IFC3, N_atoms, D)

    IFC3 = apply_tols!(IFC3,tol)

    return DenseForceConstants(IFC3, energy_unit(pot) / length_unit(pot)^3, tol)

end