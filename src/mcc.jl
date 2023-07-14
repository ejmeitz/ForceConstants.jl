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
    d = device!(device_id)
    @info "Using $d to compute MCC3\n"
    N_modes = size(phi)[1]

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
    
    #Apply tolerances -- this allocates maybe just loop?
    K3[abs.(K3) .< tol] .= 0.0
        
    return K3

end


# """
# Takes flat version of MCC3 matrix and converts it back to an (N x N x N) symmetric tensor
# using the tetrahedral numbers.
# """
# function unroll_MCC3(K3::AbstractVector, N_modes)
#     K3 = zeros(N_modes, N_modes, N_modes)
# end


# #TODO: Do I need to transpose phi to get coalesce??
# function to_mcc2(Ψ_mw_sparse::ThirdOrderSparse, phi, tol)

#     N_modes = size(phi)[1]

#     # Re-organize unqiue part of K3 into 1D
#     len_K3 = length(with_replacement_combinations(1:N_modes, 3)) # N_modes choose 3
#     K3_flat = CUDA.zeros(Float32, len_K3)
  
#     for el in Ψ_mw_sparse.values
#         @cuda @views kernel2(K3_flat, el.val, phi[el.i,:], phi[el.j,:], phi[el.k,:], N_modes, len_K3)
#     end

# end

# function to_1D_index(n,m,o, N_Modes)

# function kernel2(K3_flat, Ψ_value, phi_i, phi_j, phi_k, N_modes, maxDim)
#     n = blockIdx().x
#     m = blockIdx().y
#     o = blockIdx().z

#     idx = to_1D_index(n, m, o, N_modes)

#     if idx < maxDim #check divergence
#         K3_flat[idx] += Ψ_value * phi_i[n] * phi_j[m] * phi_k[o]
#     end
# end

# function mcc_kernel(Ψ_mw_values::Vector{F3_val})

#     threads_per_block = blockDim().x
#     idx = (blockIdx().x * threads_per_block) + threadIdx().x #one thread for every value in Ψ

#     Ψ_ijk = Ψ_mw_values[idx].val

#     shared = CuDynamicSharedArray(T, (threads_per_block,))

#     K3[n,m,l] += Ψ_ijk.val * phi[Ψ_ijk.i,n] * phi[Ψ_ijk.j,m] * phi[Ψ_ijk.k,o]


# end