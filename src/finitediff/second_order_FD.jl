export second_order, second_order_finite_diff_single

function second_order(sys_eq::SuperCellSystem{D}, pot::Potential,
    calc::FiniteDiffCalculator) where D

    N_atoms = n_atoms(sys_eq)
    IFC2 = zeros(D*N_atoms, D*N_atoms)

    r_cut_sq = calc.r_cut^2

    Threads.@threads for i in range(1, N_atoms)
        rᵢⱼ = similar(sys_eq.atoms.position[1])
        for j in range(i+1, N_atoms)

            rᵢⱼ .= sys_eq.atoms.position[i] .- sys_eq.atoms.position[j]
            nearest_mirror!(rᵢⱼ, sys_eq.box_sizes_SC)
            dist_ij_sq = sum(x -> x^2, rᵢⱼ)

            if dist_ij_sq < r_cut_sq
                for α in range(1,D)
                    for β in range(1,D)
                        ii = D*(i-1) + α
                        jj = D*(j-1) + β
                        IFC2[ii,jj] = ustrip(second_order_finite_diff_single(sys_eq, pot, [i,j], [α,β],
                            calc.r_cut, calc.h))
                        IFC2[jj,ii] = IFC2[ii,jj]
                    end
                end
            end

        end 
    end

    # Threads.@threads for i in range(1, N_atoms) # index of block matrix
    #     for α in range(1,D)
    #         for β in range(1,D)
    #             ii = D*(i-1) + α
    #             jj = D*(i-1) + β # i == j because we're on diagonal
    #             IFC2[ii,jj] = -1*sum(IFC2[ii, β:D:end])
    #         end
    #     end
    # end


    return DenseForceConstants(IFC2, energy_unit(pot) / length_unit(pot)^2, 0.0)

end

function second_order_finite_diff_single(sys_eq::SuperCellSystem{3}, pot::Potential, atom_idxs,
     cartesian_idxs, r_cut, h)

   h = uconvert(length_unit(pot), h)
   N_atoms = n_atoms(sys_eq)


#    @assert length(atom_idxs) == 2
#    @assert length(cartesian_idxs) == 2
   @assert all(atom_idxs .<= N_atoms) && all(atom_idxs .>= 1) "Atom indexes out of range, must be in 1:$(N_atoms)"
   @assert all(cartesian_idxs .<= 3) && all(cartesian_idxs .>= 1) "Cartesian indices must be 1, 2, or 3"

   posns = positions(sys_eq)


#    force_unit = energy_unit(pot) / length_unit(pot)
#    Fᵦ₀ = force_loop_j(pot, posns, force_unit, r_cut, sys_eq.box_sizes_SC, N_atoms, atom_idxs[2], cartesian_idxs[2])
#    posns[atom_idxs[1]][cartesian_idxs[1]] += h
#    ΔFᵦ = force_loop_j(pot, posns, force_unit, r_cut, sys_eq.box_sizes_SC, N_atoms, atom_idxs[2], cartesian_idxs[2])
#    posns[atom_idxs[1]][cartesian_idxs[1]] -= h
#    return - (ΔFᵦ - Fᵦ₀)/h 



    if atom_idxs[1] == atom_idxs[2] && cartesian_idxs[1] == cartesian_idxs[2]
        energies = zeros(3)*energy_unit(pot)
        combos = [h,0.0*length_unit(pot),-h]

        for (c,combo) in enumerate(combos)
            posns[atom_idxs[1]][cartesian_idxs[1]] += combo
            # posns[atom_idxs[1]][cartesian_idxs[2]] += combo

            energies[c] = energy_loop(pot, posns, sys_eq.box_sizes_SC, N_atoms, r_cut)

            posns[atom_idxs[1]][cartesian_idxs[1]] -= combo
            # posns[atom_idxs[1]][cartesian_idxs[2]] -= combo
        end

        return (1/(h^2))*(energies[1] - 2*energies[2] + energies[3])
    else
        energies = zeros(4)*energy_unit(pot)
        combos = [[h,h],[-h,-h],[h,-h],[-h,h]]

        for (c,combo) in enumerate(combos)

            posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]
            posns[atom_idxs[2]][cartesian_idxs[2]] += combo[2]

            energies[c] = energy_loop(pot, posns, sys_eq.box_sizes_SC, N_atoms, r_cut)

            #Un-modify
            posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
            posns[atom_idxs[2]][cartesian_idxs[2]] -= combo[2]
        end

        return (1/(4*h^2))*(energies[1] + energies[2] - energies[3] - energies[4])
    end
end