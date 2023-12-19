export third_order

function third_order(sys_eq::SuperCellSystem{3}, pot::Potential,
    calc::FiniteDiffCalculator)

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