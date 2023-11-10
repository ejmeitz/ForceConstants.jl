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

    IFC2 = apply_tols!(IFC2,tol)

    return SecondOrderMatrix(IFC2, unit(pot.ϵ / pot.σ^2), tol)

end

function second_order_AD_test(sys::SuperCellSystem{D}, pot::StillingerWeberSilicon, tol; r_cut = pot.r_cut) where D
    vars = make_variables(:r, D)
    r_norm = sqrt(sum(x -> x^2, vars))
    println("at start")
    #2nd order part 
    pot2_symbolic = pair_potential_nounits(pot, r_norm)
    H2_symbolic = hessian(pot2_symbolic, vars)
    H2_exec = make_function(H2_symbolic, vars)

    println("here")
    #3rd order part(s)
    pot3_symbolic = threebody_potential_nounits(pot, r_norm)
    H3_symbolic = hessian(pot3_symbolic, vars)
    H3_exec = make_function(H3_symbolic, vars)
    println("here3")
    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)

    for i in range(1, N_atoms)
        for j in range(i + 1, N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)

            if norm(rᵢⱼ) < r_cut

                #2nd order part
                ij_block = H2_exec(ustrip.(rᵢⱼ))

                IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(j-1) + 1 : D*(j-1) + D] .+= ij_block
                IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(i-1) + 1 : D*(i-1) + D] .+= ij_block

                #3rd order part(s)
                for k in range(j + 1, N_atoms)
                    rᵢₖ = sys.atoms.position[i] .- sys.atoms.position[k]
                    nearest_mirror!(rᵢₖ, sys.box_sizes_SC)
                    if norm(rᵢₖ) < r_cut
                        rⱼₖ = sys.atoms.position[j] .- sys.atoms.position[k]
                        nearest_mirror!(rⱼₖ, sys.box_sizes_SC)

                        ij_block = H3_exec(ustrip.(rᵢⱼ))
                        IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(j-1) + 1 : D*(j-1) + D] .+= ij_block
                        IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(i-1) + 1 : D*(i-1) + D] .+= ij_block

                        ik_block = H3_exec(ustrip.(rᵢₖ))
                        IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(k-1) + 1 : D*(k-1) + D] .+= ik_block
                        IFC2[D*(k-1) + 1 : D*(k-1) + D, D*(i-1) + 1 : D*(i-1) + D] .+= ik_block

                        jk_block = H3_exec(ustrip.(rⱼₖ))
                        IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(k-1) + 1 : D*(k-1) + D] .+= jk_block
                        IFC2[D*(k-1) + 1 : D*(k-1) + D, D*(j-1) + 1 : D*(j-1) + D] .+= jk_block
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

    return ThirdOrderMatrix(IFC3, unit(pot.ϵ / pot.σ^3), tol)

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