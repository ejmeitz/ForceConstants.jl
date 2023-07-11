
#could remove allocations here and just make it a ! function
function nearest_mirror(r_ij, box_sizes)
    r_x = r_ij[1]; r_y = r_ij[2]; r_z = r_ij[3]
    L_x, L_y, L_z = box_sizes
    if r_x > L_x/2
        r_x -= L_x
    elseif r_x < -L_x/2
        r_x += L_x
    end
        
    if r_y > L_y/2
        r_y -= L_y
    elseif r_y < -L_y/2
        r_y += L_y  
    end
        
    if r_z > L_z/2
        r_z -= L_z
    elseif r_z < -L_z/2
        r_z += L_z
    end
    
    return [r_x,r_y,r_z] 
end


# struct F2_val
#     i::Int32
#     j::Int32
#     val::Float32
# end

# function COO_to_ELL(Φ::Vector{F2_val})

#     IS = [];
#     JS = [];
#     vals = [];
#     for el in Φ
#         push!(IS, el.i); push!(JS, el.j)
#         push!(vals, el.val)
#     end

#     unique_rows = unique(IS)
#     max_els_per_row = maximum([count(==(i), IS) for i in unique_rows])

#     ELL_matrix = zeros(length(unique_rows), max_els_per_row)
#     ELL_column_data = similar(ELL_matrix)

#     row_counts = ones(length(IS))
#     for el in Φ
#         ELL_matrix[el.i, Int(row_counts[el.i])] = el.val
#         ELL_column_data[el.i, Int(row_counts[el.i])] = el.j
#         row_counts[el.i] += 1
#     end

#     return ELL_matrix, ELL_column_data
# end

# function to_sparse(F2::Array{T,2}) where T
#     num_nonzero = sum(F2 .!= 0.0)
#     F2_non_zero = Vector{F2_val}(undef,(num_nonzero,))
#     N = size(F2)[1]

#     count = 1
#     for i in range(1,N)
#         for j in range(1,N)
#             if F2[i,j] != 0
#                 F2_non_zero[count] = F2_val(i, j,F2[i,j])
#                 count += 1
#             end
#         end
#     end
#     return F2_non_zero
# end