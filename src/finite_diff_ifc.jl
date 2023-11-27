export check_ifc, second_order_finite_diff, third_order_finite_diff, fourth_order_finite_diff


function second_order_finite_diff(sys_eq::SuperCellSystem{3}, pot::PairPotential, atom_idxs, cartesian_idxs;
    r_cut = pot.r_cut, h = 0.04*0.5291772109u"Å")

   h = uconvert(unit(pot.σ), h)
   N_atoms = n_atoms(sys_eq)


   @assert length(atom_idxs) == 2
   @assert length(cartesian_idxs) == 2
   @assert all(atom_idxs .<= N_atoms) && all(atom_idxs .>= 1) "Atom indexes out of range, must be in 1:$(N_atoms)"
   @assert all(cartesian_idxs .<= 3) && all(cartesian_idxs .>= 1) "Cartesian indices must be 1, 2, or 3"

   posns = positions(sys_eq)

   if atom_idxs[1] != atom_idxs[2]

        energies = zeros(4)*energy_unit(pot)
        combos = [[h,h],[-h,-h],[h,-h],[-h,h]]

        for (c,combo) in enumerate(combos)

            posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]
            posns[atom_idxs[2]][cartesian_idxs[2]] += combo[2]

            energies[c] = energy_loop(pot, posns, r_cut, sys_eq.box_sizes_SC, N_atoms)

            #Un-modify
            posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
            posns[atom_idxs[2]][cartesian_idxs[2]] -= combo[2]
        end

        return (1/(4*h^2))*(energies[1] + energies[2] - energies[3] - energies[4])
    else
        energies = zeros(3)*energy_unit(pot)
        combos = [h,zero(pot.σ),-h]

        for (c,combo) in enumerate(combos)
            posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]

            energies[c] = energy_loop(pot, posns, r_cut, sys_eq.box_sizes_SC, N_atoms)

            posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
        end

        return (1/(h^2))*(energies[1] - 2*energies[2] + energies[3])

    end
end


function third_order_finite_diff(sys_eq::SuperCellSystem{3}, pot::PairPotential, atom_idxs, cartesian_idxs;
    r_cut = pot.r_cut, h = 0.04*0.5291772109u"Å")

   h = uconvert(unit(pot.σ), h)

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

#https://arxiv.org/pdf/2104.04895.pdf  Eqn 9
function fourth_order_finite_diff(sys_eq::SuperCellSystem{3}, pot::PairPotential, atom_idxs, cartesian_idxs;
    r_cut = pot.r_cut, h = 0.04*0.5291772109u"Å")

   h = uconvert(unit(pot.σ), h)

   @assert length(atom_idxs) == 4
   @assert length(cartesian_idxs) == 4
   @assert !all(atom_idxs .== atom_idxs[1]) "Cannot check self terms with finite difference"

   N_atoms = n_atoms(sys_eq)

   combos = [[h,h,h],[h,h,-h],[h,-h,h],[h,-h,-h],[-h,h,h],[-h,h,-h],[-h,-h,h],[-h,-h,-h]]
   
   force_unit = unit(force(pot,r_cut)) 
   forces = zeros(8)*force_unit

   posns = positions(sys_eq)
   
   for (c,combo) in enumerate(combos)

       posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]
       posns[atom_idxs[2]][cartesian_idxs[2]] += combo[2]
       posns[atom_idxs[3]][cartesian_idxs[3]] += combo[3]

       l = atom_idxs[4]
       force_val = zero(force(pot,1u"Å"))  
       #Calculate force on atom l
       for i in range(1,N_atoms)
           if i != l               
               r = posns[i] .- posns[l]
               nearest_mirror!(r, sys_eq.box_sizes_SC)
               dist = norm(r)
               
               #Make sure mirrored particle is in cuttoff
               if dist < r_cut
                   #Get force, potential with modified potential/ force function
                   F_mag = force(pot, dist)
                   
                   r_hat = r / dist 
                   F_ij = F_mag.*r_hat

                   #Update forces and potential
                   force_val += F_ij[cartesian_idxs[4]] #& IS THSI SIGN RIGHT??
               end
           end
       end

       forces[c] = force_val

       #Un-modify
       posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
       posns[atom_idxs[2]][cartesian_idxs[2]] -= combo[2]
       posns[atom_idxs[3]][cartesian_idxs[3]] -= combo[3]
   end

   return (1/(8*(h^3)))*(forces[1] - forces[2] - forces[3] + forces[4] - forces[5] + forces[6] + forces[7] - forces[8])
end

#& CURRENTLY DIFFERS by -1 FROM MINE
function check_ifc(sys::SuperCellSystem{D}, ifc::DenseForceConstants{O}, pot::PairPotential, n_points;
    r_cut = pot.r_cut, h = 0.04*0.5291772109u"Å") where {D,O}

    N_atoms = n_atoms(sys)
    ri = rand(1:N_atoms, (n_points,O)) #random atom indexes
    rci = rand(1:D, (n_points, O)) #random cartesian indexes

    finite_diff_funcs = [second_order_finite_diff, third_order_finite_diff, fourth_order_finite_diff]
    finite_diff_func = finite_diff_funcs[O-1]

    idxs = (ri_row, rci_row) -> [D*(ri_row[o] - 1) + rci_row[o] for o in 1:O]
    ifc_actual = @views [ifc[idxs(ri[i,:], rci[i,:])...] for i in 1:n_points]
 
    ifc_fd = zeros(n_points)*ifc.units
    for i in 1:n_points
        val = finite_diff_func(sys, pot, ri[i,:], rci[i,:]; r_cut = r_cut, h = h)
        if ustrip.(abs(val)) > ifc.tol
            ifc_fd[i] = val
        end
    end

    return ifc_fd, (ifc_actual.*ifc.units)

end

function check_ifc(sys::SuperCellSystem{D}, ifc::SparseForceConstants{O}, pot::PairPotential, n_points;
    r_cut = pot.r_cut, h = 0.04*0.5291772109u"Å") where {D,O}

    random_idxs = rand(1:length(ifc), (n_points,)) 

    finite_diff_funcs = [second_order_finite_diff, third_order_finite_diff, fourth_order_finite_diff]
    finite_diff_func = finite_diff_funcs[O-1]

    ifc_fd = zeros(n_points)*ifc.units
    ifc_actual = value.(ifc[random_idxs])
 
    for p in 1:n_points
        idxs = idx(ifc[random_idxs[p]])
        cart_idxs = ((idxs .- 1) .% D) .+ 1
        atom_idxs = ((idxs .- cart_idxs) ./ D) .- 1
        ifc_fd[p] = finite_diff_func(sys, pot, atom_idxs, cart_idxs; r_cut = r_cut, h = h)
    end

    return ifc_fd,  (ifc_actual.*ifc.units)

end