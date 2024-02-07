export third_order, third_order_finite_diff_eng

function third_order(sys_eq::SuperCellSystem{D}, pot::Potential,
    calc::FiniteDiffCalculator) where D

    throw(error("Not implemented yet"))
    
end


function third_order_finite_diff_single(sys_eq::SuperCellSystem{3}, pot::PairPotential, atom_idxs, cartesian_idxs,
    r_cut, h)

   h = uconvert(length_unit(pot), h)

   @assert length(atom_idxs) == 3
   @assert length(cartesian_idxs) == 3

   N_atoms = n_atoms(sys_eq)

   z = zero(pot.σ)
   combos = [[z,h,-h],[z,h,h],[z,-h,h],[z,-h,-h]]
   
   force_unit = unit(force(pot,r_cut)) 
   forces = zeros(4)*force_unit

   posns = positions(sys_eq)

   for (c,combo) in enumerate(combos)

       posns[atom_idxs[2]][cartesian_idxs[2]] += combo[2]
       posns[atom_idxs[3]][cartesian_idxs[3]] += combo[3]

       forces[c] = force_loop_j(pot, posns, force_unit, r_cut, sys_eq.box_sizes_SC, N_atoms, atom_idxs[1], cartesian_idxs[1])

       posns[atom_idxs[2]][cartesian_idxs[2]] -= combo[2]
       posns[atom_idxs[3]][cartesian_idxs[3]] -= combo[3]
   end

   return (1/(4*(h^2)))*(forces[1] - forces[2] + forces[3] - forces[4])
end



function third_order_finite_diff_eng(sys_eq::SuperCellSystem{3}, pot::Potential, atom_idxs, cartesian_idxs;
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