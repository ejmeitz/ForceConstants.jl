# export fourth_order

# function fourth_order(sys_eq::SuperCellSystem{3}, pot::Potential,
#     calc::FiniteDiffCalculator)

#     throw(error("Not implemented yet"))
    
# end

# #https://arxiv.org/pdf/2104.04895.pdf  Eqn 9
# function fourth_order_finite_diff_single(sys_eq::SuperCellSystem{3}, pot::PairPotential, atom_idxs, cartesian_idxs,
#     r_cut, h)

#    h = uconvert(length_unit(pot), h)

#    @assert length(atom_idxs) == 4
#    @assert length(cartesian_idxs) == 4
#    @assert !all(atom_idxs .== atom_idxs[1]) "Cannot check self terms with finite difference"

#    N_atoms = n_atoms(sys_eq)

#    combos = [[h,h,h],[h,h,-h],[h,-h,h],[h,-h,-h],[-h,h,h],[-h,h,-h],[-h,-h,h],[-h,-h,-h]]
   
#    force_unit = unit(force(pot,r_cut)) 
#    forces = zeros(8)*force_unit

#    posns = positions(sys_eq)
   
#    for (c,combo) in enumerate(combos)

#        posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]
#        posns[atom_idxs[2]][cartesian_idxs[2]] += combo[2]
#        posns[atom_idxs[3]][cartesian_idxs[3]] += combo[3]

#        l = atom_idxs[4]
#        force_val = zero(force(pot,1u"Å"))  
#        #Calculate force on atom l
#        for i in range(1,N_atoms)
#            if i != l               
#                r = posns[i] .- posns[l]
#                nearest_mirror!(r, sys_eq.box_sizes_SC)
#                dist = norm(r)
               
#                #Make sure mirrored particle is in cuttoff
#                if dist < r_cut
#                    #Get force, potential with modified potential/ force function
#                    F_mag = force(pot, dist)
                   
#                    r_hat = r / dist 
#                    F_ij = F_mag.*r_hat

#                    #Update forces and potential
#                    force_val += F_ij[cartesian_idxs[4]] #& IS THSI SIGN RIGHT??
#                end
#            end
#        end

#        forces[c] = force_val

#        #Un-modify
#        posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
#        posns[atom_idxs[2]][cartesian_idxs[2]] -= combo[2]
#        posns[atom_idxs[3]][cartesian_idxs[3]] -= combo[3]
#    end

#    return (1/(8*(h^3)))*(forces[1] - forces[2] - forces[3] + forces[4] - forces[5] + forces[6] + forces[7] - forces[8])
# end
