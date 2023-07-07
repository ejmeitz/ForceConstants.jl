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
function to_mcc(Ψ_mw::ThirdOrderSparse, phi, tol, devices::Vector{CuDevice})

    N_modes = size(Ψ_mw)[1]

    #Move F3 & phi to GPU
    cuF3_sparse = CuArray(F3_sparse)
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