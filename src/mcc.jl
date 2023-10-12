export mcc3, mcc3_custom_kernel


"""
Converts third order forces constants, `Ψ` into third order modal coupling constants (MCC).
Does not divide MCC calculation into smaller chunks. This might exhaust GPU memory.
"""
function mcc3(Ψ::CuArray{Float32, 3}, phi::CuArray{Float32, 2}, tol::Float64)

    K3 = CUDA.zeros(Float32, size(Ψ));

    @tensor begin
        K3[n,m,l] = Ψ[i,j,k] * phi[i,n] * phi[j,m] * phi[k,l]
    end

    K3_CPU = Array(K3)
    return apply_tols!(K3_CPU, tol)

end

"""
Converts third order forces constants, `Ψ` into third order modal coupling constants (MCC). The
parameter `block_size` specifies problem size when calculating the MCC. For example, if 
`block_size` is 100, the block MCC will be calculated in 100x100x100 blocks to save GPU memory.
The `phi` matrix will be automatically truncated to adjust for this.
    Try to maximize `block_size` and make it a power of 2.
"""
function mcc3(Ψ::CuArray{Float32, 3}, phi::CuArray{Float32, 2}, block_size::Int, tol::Float64)

    @assert size(phi)[1] % block_size == 0
    @assert block_size > 0

    n_blocks_per_dim = Int(size(phi)[1] / block_size)

    #Keep large storage on CPU
    K3_CPU = zeros(Float32, size(Ψ))
    tmp = zeros(Float32, (block_size,block_size,block_size))
    #Pre-allocate GPU storage to re-use
    K3_GPU_block = CUDA.zeros(Float32, (block_size, block_size, block_size))

    for i in 1:n_blocks_per_dim
        for j in 1:n_blocks_per_dim
            for k in 1:n_blocks_per_dim
                dim1_range = (block_size*(i-1) + 1):(block_size*i)
                dim2_range = (block_size*(j-1) + 1):(block_size*j)
                dim3_range = (block_size*(k-1) + 1):(block_size*k)
                K3_GPU_block =  mcc3_blocked!(K3_GPU_block, Ψ, phi, [dim1_range, dim2_range, dim3_range])
                
                copyto!(tmp, K3_GPU_block) #TODO can i remove this? uses scalar indexing to copy directly into big array
                copyto!(view(K3_CPU, dim1_range, dim2_range, dim3_range), tmp)
            end
        end
    end

    #TODO this calculates more than it needs to atm, can only calculate unique part and rotate 
    #TODO to enforce symmetry, kinda annoying to implement, will matter more when system > GPU RAM

    return apply_tols!(K3_CPU, tol)
end


function mcc3_blocked!(K3::CuArray{Float32, 3}, Ψ::CuArray{Float32, 3}, phi::CuArray{Float32, 2}, idx_ranges::Vector{UnitRange{Int64}})

    phi_block1 = view(phi, :, idx_ranges[1])
    phi_block2 = view(phi, :, idx_ranges[2])
    phi_block3 = view(phi, :, idx_ranges[3])

    @tensor begin
        K3[n,m,l] = Ψ[i,j,k] * phi_block1[i,n] * phi_block2[j,m] * phi_block3[k,l]
    end

    return K3
end


const BLOCK_SIZE = 128

#Each thread calculates fiber of K_nml
function mcc3_custom_kernel(num_psi_idxs::Int32, mcc::CuDeviceArray{Float32,1}, phi::CuDeviceArray{Float32,2},
         Ψ_vals::CuDeviceArray, is::CuDeviceArray{Int32,1}, js::CuDeviceArray{Int32,1}, ks::CuDeviceArray{Int32,1}) 

    n = (blockIdx().x-1i32) * blockDim().x + threadIdx().x
    m = (blockIdx().y-1i32) * blockDim().y + threadIdx().y
    linearThreadIdx = (blockDim().x * threadIdx().y) + threadIdx().x #& OFF BY 1 iSSUES?

    #To store rows of phi matrix needed by current values of Psi
    phi_i = CuStaticSharedArray(Float32, BLOCK_SIZE)
    phi_j = CuStaticSharedArray(Float32, BLOCK_SIZE)
    phi_k = CuStaticSharedArray(Float32, BLOCK_SIZE)


    # Treat thread index as index into Ψ now and accumulate values into MCC
    for i in 1:BLOCK_SIZE:num_psi_idxs
        psi_idx = i + linearThreadIdx

        if psi_idx <= num_psi_idxs
            #Move Phi Data Into Shared Memory
            phi_i[linearThreadIdx] = phi[is[psi_idx]]
            phi_j[linearThreadIdx] = phi[js[psi_idx]]
            phi_k[linearThreadIdx] = phi[ks[psi_idx]] #& This might be the only one that needs to be in shared mem
            sync_threads()

            #Calculate fiber of MCC
            for l in 1:N_modes
                mcc[n,m,l] += Ψ_vals[psi_idx]*phi_i[n]*phi_j[m]*phi_k[l]
            end
            sync_threads()

        end
    end

end

# const BLOCK_SIDE = 32
# const GRID_SIDE = ceil(sqrt(N)/BLOCK_SIDE)
# @cuda blocks = (GRID_SIDE,GRID_SIDE) threads = (BLOCK_SIDE,BLOCK_SIDE) mcc3_custom_kernel(

# )