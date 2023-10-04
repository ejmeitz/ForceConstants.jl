export mcc3


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

