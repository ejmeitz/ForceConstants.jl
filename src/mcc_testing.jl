

@inbounds  function mcc3_kernel!(K3, rows, cols, Ψ_ijk, phi_i, phi_j, phi_k, maxIdx)
    thread_id = ((blockIdx().x - 1) * blockDim().x) + threadIdx().x;  

    if thread_id <= maxIdx
        idx1 = rows[thread_id]
        idx2 = cols[thread_id]
        idx3 =  row_col_2_depth(idx1, idx2, thread_id)
        K3[thread_id] += Ψ_ijk * phi_i[idx1] * phi_j[idx2] * phi_k[idx3]
    end

    return
end

const CUDA_MAX_THREADS_PER_BLOCK = CUDA.attribute(device(), CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK)
const MAX_BLOCKS_PER_DIM_X = CUDA.attribute(device(), CUDA.DEVICE_ATTRIBUTE_MAX_GRID_DIM_X)
const WARP_SIZE = CUDA.attribute(device(),CUDA.DEVICE_ATTRIBUTE_WARP_SIZE)
function mcc3(Ψ, cuPhi, N_modes)

    M = tetra_num(N_modes)
    @assert M < typemax(Int32) "Number of elements in MCC3 overflows Int32"
    K3 = CUDA.zeros(M)
    THREADS_PER_BLOCK, NUM_BLOCKS = calculate_CUDA_dimension(N_modes, M)

    # How to launch thread blocks such that phi is accessed in order??
    # need shenanegins with warp stuff and shared memory

    #Do I index Ψ inside of kernel? will this help get better memory access patterns??

    #Is the overhead from launching kerenls important??
    
    for Ψ_ijk in Ψ
        phi_i = view(cuPhi, Ψ_ijk.i, :) #this accesses in not-col major #TODO
        phi_j = view(cuPhi, Ψ_ijk.j, :)
        phi_k = view(cuPhi, Ψ_ijk.k, :)
        @cuda threads=THREADS_PER_BLOCK blocks=NUM_BLOCKS mcc3_kernel!(K3, rows, cols, Ψ_ijk, phi_i, phi_j, phi_k, M)
    end

end

function calculate_CUDA_dimension(N_modes, M)
    THREADS_PER_BLOCK = min(N_modes, CUDA_MAX_THREADS_PER_BLOCK)
    NUM_BLOCKS = cld(M, THREADS_PER_BLOCK)

    @assert NUM_BLOCKS < MAX_BLOCKS_PER_DIM_X
    @assert (NUM_BLOCKS*THREADS_PER_BLOCK) >= M

    return THREADS_PER_BLOCK, NUM_BLOCKS
end