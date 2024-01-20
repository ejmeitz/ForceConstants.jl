# export check_ifc

# function check_ifc(sys::SuperCellSystem{D}, ifc::DenseForceConstants{O}, pot::PairPotential, n_points;
#     calc::FiniteDiffCalculator) where {D,O}

#     if calc.tol != ifc.tol
#         @warn "IFC was calculated with different tolerance than specified in calculator. 
#                 Using calculator tolerance."
#     end

#     N_atoms = n_atoms(sys)
#     ri = rand(1:N_atoms, (n_points,O)) #random atom indexes
#     rci = rand(1:D, (n_points, O)) #random cartesian indexes

#     finite_diff_funcs = [second_order_finite_diff_single, third_order_finite_diff_single, fourth_order_finite_diff_single]
#     finite_diff_func = finite_diff_funcs[O-1]

#     idxs = (ri_row, rci_row) -> [D*(ri_row[o] - 1) + rci_row[o] for o in 1:O]
#     ifc_actual = @views [ifc[idxs(ri[i,:], rci[i,:])...] for i in 1:n_points]
 
#     ifc_fd = zeros(n_points)*ifc.units
#     for i in 1:n_points
#         val = finite_diff_func(sys, pot, ri[i,:], rci[i,:]; calc.r_cut, calc.h)
#         if ustrip.(abs(val)) > ifc.tol
#             ifc_fd[i] = val
#         end
#     end

#     return ifc_fd, (ifc_actual.*ifc.units)

# end

# function check_ifc(sys::SuperCellSystem{D}, ifc::SparseForceConstants{O}, pot::PairPotential, n_points;
#     calc::FiniteDiffCalculator) where {D,O}

#     if calc.tol != ifc.tol
#         @warn "IFC was calculated with different tolerance than specified in calculator. 
#                 Using calculator tolerance."
#     end

#     random_idxs = rand(1:length(ifc), (n_points,)) 

#     finite_diff_funcs = [second_order_finite_diff_single, third_order_finite_diff_single, fourth_order_finite_diff_single]
#     finite_diff_func = finite_diff_funcs[O-1]

#     ifc_fd = zeros(n_points)*ifc.units
#     ifc_actual = value.(ifc[random_idxs])
 
#     for p in 1:n_points
#         idxs = idx(ifc[random_idxs[p]])
#         cart_idxs = ((idxs .- 1) .% D) .+ 1
#         atom_idxs = ((idxs .- cart_idxs) ./ D) .- 1
#         ifc_fd[p] = finite_diff_func(sys, pot, atom_idxs, cart_idxs, calc.r_cut, calc.h)
#     end

#     return ifc_fd,  (ifc_actual.*ifc.units)

# end
