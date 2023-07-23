

@inbounds  function mcc3_kernel!(K3, rows, cols, Ψ_ijk, phi_i, phi_j, phi_k)
    i = ((blockIdx().x - 1) * blockDim().x) + threadIdx().x;

    idx1 = rows[i]
    idx2 = cols[i]
    idx3 =  row_col_2_depth(idx1, idx2, i) 

    K3[i] += Ψ_ijk * phi_i[idx1] * phi_j[idx2] * phi_k[idx3]

    return
end

const CUDA_MAX_THREADS_PER_BLOCK = 1024
const WARP_SIZE = 32
function mcc3(Ψ, phi, N_modes)

    THREADS_PER_BLOCK = min(N_modes, CUDA_MAX_THREADS_PER_BLOCK)
    M = tetra_num(N_modes)
    K3 = zeros(Float32, M)

    # How to launch thread blocks such that phi is accessed in order??
    # need shenanegins with warp stuff and shared memory
    for Ψ_ijk in Ψ
        phi_i = view(phi, Ψ_ijk.i, :) #this accesses in not-col major #TODO
        phi_j = view(phi, Ψ_ijk.j, :)
        phi_k = view(phi, Ψ_ijk.k, :)
        @cuda threads=THREADS_PER_BLOCK blocks=cld(M,THREADS_PER_BLOCK) mcc3_kernel!(K3, rows, cols, Ψ_ijk, phi_i, phi_j, phi_k)
    end

end