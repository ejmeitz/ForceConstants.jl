export to_mcc

"""
This file constains functions for converting force constants
into modal coupling constants. These functions are extremely expensive
beyond second order and require GPU resources to be efficient.
"""

### Second Order Functions ###
function to_mcc(dynmat::SecondOrderMatrix, num_rigid_translation = 3)
    freqs_sq, _ = get_modes(dynmat, num_rigid_translation)
    return freqs_sq
end

### Third Order Functions ###

#TODO: Custom Kernel, kernel maybe should be all of the operations
function gpu_k3_kernel(cuF3_sparse_mw, cuPhi1, cuPhi2, cuPhi3)
    f = (f3_data) -> f3_data.val * cuPhi1[f3_data.i] * cuPhi2[f3_data.j] * cuPhi3[f3_data.k]
   return mapreduce(f, +, cuF3_sparse_mw)
end

"""
Takes in sparse third order data that has already been mass weighted

Parameters:
 - Ψ_mw: Sparse, mass-weighted third order force constants
 - phi: Eigenvectors (mode shapes) for the system
 - tol: Values less than tol will be set to 0
 - devices: List of GPUs to target
"""
function to_mcc(Ψ_mw_sparse::ThirdOrderSparse, phi, tol, device_id)#; devices::Vector{CuDevice} = nothing)
    device!(device_id)
    @info "Using $(device_id) to compute MCC3\n"
    N_modes = length(Ψ_mw_sparse)

    #Move F3 & phi to GPU
    cuF3_sparse = CuArray(Ψ_mw_sparse.values)
    cuPhi = CuArray(phi)


    
    K3 = zeros(N_modes,N_modes,N_modes)
    K3_indices = with_replacement_combinations(range(1,N_modes), 3)
    
    for idx in K3_indices
        m,n,o = idx
        
        val = @views gpu_k3_kernel(cuF3_sparse, cuPhi[:,m], cuPhi[:,n], cuPhi[:,o])
        K3[m,n,o] = val
        K3[m,o,n] = val
        K3[o,n,m] = val
        K3[o,m,n] = val
        K3[n,m,o] = val
        K3[n,o,m] = val
    end
    
    #Apply tolerances
    K3[abs.(K3) .< tol] .= 0.0
        
    return K3

end

# function to_mcc(Ψ_mw::ThirdOrderSparse)

#     # (I, J, K, V) = findnz(B)

#     #TODO: is it faster to just have arrays of indices vs using F3_val struct??
#     #NOTE THIS ONLY DOES UPPER TRIANGLE OF K3

#     K3_indices = with_replacement_combinations(range(1, N_modes), 3)
#     for i in eachindex(Ψ_mw.values)
#         begin
#             for o in 1:N_modes
#                 for n in 1:o
#                     for m in 1:n
#                         K3[m, n, o] += Ψ_mw[i].val * phi[Ψ_mw[i].i, m] * phi[Ψ_mw[i].j, n] * phi[Ψ_mw[i].k, o]
#                     end
#                 end
#             end
#         end
#         synchronize()
#     end

#     for o in 1:N_modes
#         for n in 1:o
#             for m in 1:n
#                 K3[m, n, o] = sum([(val * phi[i,m] * phi[j,n] * phi[k,o]) for _ in eachindex(Ψ_mw.values)]) 
#             end
#         end
#     end



# end

# function to_mcc2(Ψ_mw_sparse::ThirdOrderSparse, phi, tol)

#     for el in Ψ_mw_sparse.values
#         @cuda kernel2()
#     end

# end

# function kernel2()


# end

# function mcc_kernel(Ψ_mw_values::Vector{F3_val})

#     threads_per_block = blockDim().x
#     idx = (blockIdx().x * threads_per_block) + threadIdx().x #one thread for every value in Ψ

#     Ψ_ijk = Ψ_mw_values[idx].val

#     shared = CuDynamicSharedArray(T, (threads_per_block,))

#     K3[n,m,l] += Ψ_ijk.val * phi[Ψ_ijk.i,n] * phi[Ψ_ijk.j,m] * phi[Ψ_ijk.k,o]


# end