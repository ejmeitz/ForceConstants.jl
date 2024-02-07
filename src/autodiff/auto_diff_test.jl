export calculate_force_constants, calc_derivs2, calc_derivs3, make_test_sys,third_order_finite_diff_test

function calc_derivs2(pot, D)

    return three_body_second_derivs(pot, D)

end

function calc_derivs3(pot, D)

    return three_body_third_derivs(pot, D)

end

function make_test_sys()

    r = [[0.0, 0.0, 0.0],
        [2.0,1.0,1.0],
        [1.0,2.0,1.0]]

    #Put in center of box away from edges
    r[1] .+= 5.0
    r[2] .+= 5.0
    r[3] .+= 5.0

    masses = [1,1,1]

    return SuperCellSystem(r*u"Å", masses, [20,20,20]*u"Å")

end

function calculate_force_constants(sys, derivs2, derivs3)

    r = positions(sys)

    dist_12 = norm(r[1] .- r[2])
    dist_13 = norm(r[1] .- r[3])
    dist_23 = norm(r[2] .- r[3])

    pot = StillingerWeberSilicon()

    @assert all([dist_12, dist_13, dist_23] .< pot.r_cut)
    D = 3

    H2_exec, H2_exec_ij, H2_exec_ik = derivs2
        #three_body_second_derivs(pot, D)

    H3_exec_ij, H3_exec_iij, H3_exec_iik, H3_exec_ijj, H3_exec_ijk, H3_exec_ikk = derivs3
        #three_body_third_derivs(pot, D)

    IFC2 = zeros(9,9)
    IFC3 = zeros(9,9,9)

    block2 = zeros(3,3)
    block3 = zeros(3,3,3)
    for i in range(1,3)
        i_rng = D*(i-1) + 1 : D*(i-1) + D
        for j in range(1,3)
            j_rng = D*(j-1) + 1 : D*(j-1) + D
            if i != j
                r_ij = r[i] .- r[j]

                if i < j
                    block2 .= -H2_exec(ustrip.(r_ij))
                    IFC2[i_rng, j_rng] .+= block2
                    IFC2[j_rng, i_rng] .+= block2

                    block3 .= -H3_exec_ij(ustrip.(r_ij))
                    set_third_order_terms!(IFC3, i_rng, j_rng, block3)
                    set_third_order_terms!(IFC3, j_rng, i_rng, -block3)
                end

                for k in range(j+1,3)
                    if i != k
                        r_ik = r[i] .- r[k]
                        r_j = r[i] .- r_ij
                        r_k = r[i] .- r_ik
                        k_rng = D*(k-1) + 1 : D*(k-1) + D

                        r_arr = ustrip.([r[i]; r_j; r_k])

                        block2 .= H2_exec_ij(r_arr)
                        IFC2[i_rng, j_rng] .+= block2
                        IFC2[j_rng, i_rng] .+= block2
                        
                        block2 .= H2_exec_ik(r_arr)
                        IFC2[i_rng, k_rng] .+= block2
                        IFC2[k_rng, i_rng] .+= block2

                        block3 .= H3_exec_iij(r_arr)
                        set_third_order_terms!(IFC3, i_rng, j_rng, block3)

                        block3 .= H3_exec_ijj(r_arr)
                        set_third_order_terms!(IFC3, j_rng, i_rng, block3)

                        block3 .= H3_exec_iik(r_arr)
                        set_third_order_terms!(IFC3, i_rng, k_rng, block3)

                        block3 .= H3_exec_ikk(r_arr)
                        set_third_order_terms!(IFC3, k_rng, i_rng, block3)

                        block3 .= H3_exec_ijk(r_arr)
                        set_third_order_terms!(IFC3, i_rng, j_rng, k_rng, block3)
                    end
                end
            end

        end
    end

    return IFC2, IFC3
end



function third_order_finite_diff_test(sys_eq::SuperCellSystem{3}, pot::Potential, atom_idxs, cartesian_idxs;
    r_cut = pot.r_cut, h = 0.04*0.5291772109u"Å")

    h = uconvert(length_unit(pot), h)
    N_atoms = n_atoms(sys_eq)


    @assert length(atom_idxs) == 3
    @assert length(cartesian_idxs) == 3
    @assert all(atom_idxs .<= N_atoms) && all(atom_idxs .>= 1) "Atom indexes out of range, must be in 1:$(N_atoms)"
    @assert all(cartesian_idxs .<= 3) && all(cartesian_idxs .>= 1) "Cartesian indices must be 1, 2, or 3"

    posns = [Vector(a) for a in positions(sys_eq)]
    z = zero(0.0*length_unit(pot))

    if (atom_idxs[1] != atom_idxs[2]) && (atom_idxs[2] != atom_idxs[3]) && (atom_idxs[1] != atom_idxs[3])
        energies = zeros(8)*energy_unit(pot)
        combos = [[h,h,h],[h,-h,-h],[-h,-h,h],[-h,h,-h],[-h,-h,-h],[-h,h,h],[h,-h,h],[h,h,-h]]

        for (c,combo) in enumerate(combos)
            posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]
            posns[atom_idxs[2]][cartesian_idxs[2]] += combo[2]
            posns[atom_idxs[3]][cartesian_idxs[3]] += combo[3]

            energies[c] = energy_loop(pot, posns, sys_eq.box_sizes_SC, N_atoms, r_cut)

            posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
            posns[atom_idxs[2]][cartesian_idxs[2]] -= combo[2]
            posns[atom_idxs[3]][cartesian_idxs[3]] -= combo[3]
        end

        return (1/(8*(h^3)))*(energies[1] + energies[2] + energies[3] + energies[4] -
                 energies[5] - energies[6] - energies[7] - energies[8])
        
    elseif (atom_idxs[1] == atom_idxs[2]) && (atom_idxs[2] == atom_idxs[3]) && (atom_idxs[1] == atom_idxs[3])
        energies = zeros(4)*energy_unit(pot)
        combos = [[-2*h,z,z],[-h,z,z],[h,z,z],[2*h,z,z]]

        for (c,combo) in enumerate(combos)
            posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]

            energies[c] = energy_loop(pot, posns, sys_eq.box_sizes_SC, N_atoms, r_cut)

            posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
        end

        return (1/(2*(h^3)))*(-energies[1] + 2*energies[2] - 2*energies[3] + energies[4])
    elseif (atom_idxs[1] != atom_idxs[2]) && (atom_idxs[1] == atom_idxs[3])

        #* ADD PERMUTATIONS OF THIS
        energies = zeros(6)*energy_unit(pot)
        combos = [[h,h,z],[-h,h,z],[z,-h,z],[-h,-h,z],[h,-h,z],[z,h,z]]

        for (c,combo) in enumerate(combos)

            posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]
            posns[atom_idxs[2]][cartesian_idxs[2]] += combo[2]

            energies[c] = energy_loop(pot, posns, sys_eq.box_sizes_SC, N_atoms, r_cut)

            posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
            posns[atom_idxs[2]][cartesian_idxs[2]] -= combo[2]
        end

        return (1/(2*(h^3)))*(energies[1] + energies[2] + 2*energies[3] - energies[4] - energies[5] - 2*energies[6])

    else
        @warn "Shouldnt be here"
    end 
end