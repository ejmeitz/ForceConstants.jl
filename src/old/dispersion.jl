# export get_dispersion_points

# function get_dispersion_points(sys::UnitCellSystem{3}, pot::PairPotential;
#                  directions = ([1.0, 0.0, 0.0],), group_by = :branch, tol = 1e-12, unit_system = :REAL)
    
#     @assert (group_by ∈ [:branch, :k_point]) "group_by can only be :branch or :k_point"

#     Δk = (2*pi)./sys.box_sizes_SC
#     left_zone_edges = (-pi./sys.box_sizes_UC) .+ Δk

#     @assert all(sys.num_unit_cells[1] .== sys.num_unit_cells) "Current implementation assumes same number of units cells in each direction"
#     out = Dict()
#     for direction in directions
#         k_points = []
#         k_current = left_zone_edges .* direction
#         for i in range(1,sys.num_unit_cells[1]) #this is where cubic assumption matters
#             push!(k_points, copy(k_current))
#             k_current .+= (Δk .* direction)
#         end

#         ω_all = Dict()
#         for k_point in k_points
#             dynmat_uc = dynamicalMatrix(sys, pot, SVector{3}(k_point), tol)
#             ω_sq , _ = get_modes(dynmat_uc, 3)
#             if group_by == :k_point
#                 ω_all[ustrip(k_point)] = ω_sq
#             elseif group_by == :branch
#                 throw(error("Unimplemented"))
#             end
#         end
#         out[tuple(direction...)] = ω_all
#     end




#     #Sqrt and convert units on freqs
#     for b in keys(out)
#         ω_all = out[b]
#         for k in keys(ω_all)

#             if ustrip.(k) == [0.0, 0.0, 0.0] #if at origin set rigid translation modes to 0
#                 idxs = sortperm(abs.(ω_all[k]))
#                 ω_all[k][idxs[1:3]] .= 0.0
#             end

#             if unit_system == :REAL
#                 ω_all[k] = sqrt.(ω_all[k] .* 4184 .* 1000 .* 1e20)./(2*pi*1e12)
#             else
#                 throw(error("Only real units supported at this time"))
#             end
#         end
#     end

#     return out
# end