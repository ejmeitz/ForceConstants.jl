export energy_loop, force_loop

#Kronicker Delta
δₖ(x,y) = ==(x,y)


function nearest_mirror!(r_ij, box_sizes)
    for i in eachindex(box_sizes)  
        if r_ij[i] > box_sizes[i]/2
            r_ij[i] -= box_sizes[i]
        elseif r_ij[i] < -box_sizes[i]/2
            r_ij[i] += box_sizes[i]
        end
    end
    
    return r_ij
end


function apply_tols!(arr, tol)
    Threads.@threads for i in eachindex(arr)
        if abs(arr[i]) < tol
            arr[i] = 0.0
        end
    end
    return arr
end

function energy_loop(pot::PairPotential, posns, eng_unit, r_cut, box_sizes, N_atoms)

    U_total = 0.0*eng_unit

    for i in range(1,N_atoms)
        for j in range(i+1, N_atoms)             
            r = posns[i] .- posns[j]
            nearest_mirror!(r, box_sizes)
            dist = norm(r)
            if dist < r_cut
                U_total += potential(pot,dist)
            end
        end
    end

    return U_total

end

function energy_loop(pot::StillingerWeberSilicon, posns, eng_unit, r_cut, box_sizes, N_atoms)

    U_total = 0.0*eng_unit

    for i in range(1,N_atoms)
        for j in range(1, N_atoms)

            rᵢⱼ = posns[i] .- posns[j]
            nearest_mirror!(rᵢⱼ, box_sizes)
            dist_ij = norm(rᵢⱼ) 

            if i < j && dist_ij < r_cut
                U_total += pair_potential(pot, rᵢⱼ)
            end
            
            for k in range(j+1, N_atoms)  
                # i central
                rᵢₖ = posns[i] .- posns[k]
                nearest_mirror!(rᵢₖ, box_sizes)
                dist_ik = norm(rᵢₖ)    
                
                # if dist_ij < r_cut && dist_ik < r_cut 
                U_total += three_body_potential(pot, rᵢⱼ, rᵢₖ)
                # end
            end
        end
    end
    return U_total
end

# function energy_loop(pot::StillingerWeberSilicon, posns, eng_unit, r_cut, box_sizes, N_atoms)

#     U_total = 0.0*eng_unit

#     for i in range(1,N_atoms)
#         for j in range(i+1, N_atoms)
#             rᵢⱼ = posns[i] .- posns[j]
#             nearest_mirror!(rᵢⱼ, box_sizes)
#             dist_ij = norm(rᵢⱼ) 
#             if dist_ij < r_cut
#                 U_total += pair_potential(pot, rᵢⱼ)
            
#                 for k in range(j+1, N_atoms)

#                     #i as central atom
#                     rᵢₖ = posns[i] .- posns[k]
#                     nearest_mirror!(rᵢₖ, box_sizes)
#                     dist_ik = norm(rᵢₖ)    

#                     if dist_ik < r_cut 
#                         U_total += three_body_potential(pot, rᵢⱼ, rᵢₖ)
#                     end

#                     #j,k as central atom
#                     rⱼᵢ = -rᵢⱼ
#                     rⱼₖ = rᵢₖ .- rᵢⱼ #rᵢⱼ .- rᵢₖ
#                     rₖᵢ = -rᵢₖ
#                     rₖⱼ = -rⱼₖ
#                     dist_jk = norm(rⱼₖ)
#                     if dist_jk < r_cut
#                         U_total += three_body_potential(pot, rⱼᵢ, rⱼₖ) #j central atom
#                         if dist_ik < r_cut
#                             U_total += three_body_potential(pot, rₖᵢ, rₖⱼ) #k central atom 
#                         end
#                     end
                    
#                 end
#             end
#         end
#     end
#     return U_total
# end

# Calculates force on atom j in β direction
function force_loop_j(pot::PairPotential, posns, force_unit, r_cut, box_sizes, N_atoms, j, β)

    D = length(box_sizes)
    Fⱼᵦ = zeros(D)*force_unit

    for i in range(1,N_atoms)
        if i != j            
            r = posns[i] .- posns[j]
            nearest_mirror!(r, box_sizes)
            dist = norm(r)
            if dist < r_cut
                F_mag = force(pot, dist)
                r_hat = r / dist 
                F_ij = F_mag*r_hat

                #Update forces and potential
                Fⱼᵦ += F_ij[β]
            end
        end
    end

    return Fⱼᵦ

end




# """
# nᵗʰ triangular number
# """
# function tri_num(n::Int32)
#     return Int32(0.5*n*(n+1))
# end

# """
# nᵗʰ tetrahedral number
# """
# function tetra_num(n::Int32)
#     return Int32((n*(n+1)*(n+2))/6)
# end

# """
# Takes lower part of triangular matrix and flattens it into a 1D vector.
# Returns the row of each element in the original matrix `A` and the
# flattened vector of values.
# """
# function flatten_symmetric(A::Matrix{T}) where T
#     # @assert issymmetric(A) "A must be symmetric"
#     M,N = size(A)
#     @assert M == N "All matrix dimensions must have equal length"
#     rows = zeros(Int32, tri_num(N))
#     flattened = zeros(T, tri_num(N))
#     idx = 1
#     for row in 1:M
#         for col in 1:row
#             rows[idx] = row
#             flattened[idx] = A[row, col]
#             idx += 1
#         end
#     end

#     @assert (idx-1) == length(flattened)

#     return rows, flattened
# end

# """
# Takes lower part of 3D tensor and flattens it into a 1D vector.
# Returns the row & col of each element in the original matrix `A` and the
# flattened vector of values. Assumes rows index first element, cols second element.
# """
# function flatten_symmetric(A::Array{T,3}) where T
#     M,N,O = size(A)
#     @assert ((M == N) && (N == O)) "All matrix dimensions must have equal length"
#     rows = zeros(Int32, tetra_num(N))
#     cols = zeros(Int32, tetra_num(N))
#     flattened = zeros(T, tetra_num(N))
#     idx = 1
#     for x in 1:M
#         for y in 1:x
#             for z in 1:y
#                 rows[idx] = x
#                 cols[idx] = y
#                 flattened[idx] = A[x, y, z]
#                 idx += 1
#             end
#         end
#     end

#     @assert (idx-1) == length(flattened)

#     return rows, cols, flattened
# end

# """
# Takes 1D index into flattened array and the first index (row) of the element
# in the un-flattened matrix. Returns the second index (col) of the element
# in the un-flattened matrix.
# """
# function row_2_col(idx1, idx_1D)
#     return idx_1D - tri_num(idx1 - 1)
# end

# """
# Takes 1D index into flattened array and the first and second indices of the element
# in the un-flattened tensor. Returns the third index of the element
# in the un-flattened tensor.
# """
# function row_col_2_depth(idx1, idx2, idx_1D)
#     return (idx_1D - tetra_num(idx1 - 1)) - tri_num(idx2 - 1)
# end


# function test_coordinate_generation(sz)

#     ## 2D
#     A = rand(sz,sz)
#     idx_flat = 1:tri_num(sz)
#     rows, A_flat = flatten_symmetric(A)
#     cols = row_2_col.(rows, 1:tri_num(sz))

#     for idx in idx_flat
#         @assert A_flat[idx] == A[rows[idx], cols[idx]]
#     end

#     ## 3D
#     B = rand(sz,sz,sz)
#     idx_flat = 1:tetra_num(sz)
#     rows, cols, B_flat = flatten_symmetric(B)
#     idx3 = row_col_2_depth.(rows, cols, 1:tetra_num(sz))

#     for idx in idx_flat
#         @assert B_flat[idx] == B[rows[idx], cols[idx], idx3[idx]]
#     end

# end

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