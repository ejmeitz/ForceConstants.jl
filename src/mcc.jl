export mcc3

"""
Converts third order forces constants, `Ψ` into third order modal coupling constants (MCC). The
parameter `block_size` specifies problem size when calculating the MCC. For example, if 
`block_size` is 100, the block MCC will be calculated in 100x100x100 blocks to save GPU memory.
The `phi` matrix will be automatically truncated to adjust for this. Try to minimize block_size.
"""
function mcc3(Ψ::CuArray{Float32, 3}, phi::CuArray{Float32, 2}, block_size::Int; tol = 1e-12, gpu_id::Int)
    device!(gpu_id)

    @assert size(phi)[1] % block_size == 0
    @assert block_size > 0

    n_blocks_per_dim = Int(size(phi)[1] / block_size)

    #Keep large storage on CPU
    K3 = zeros(Float32,size(Ψ))
    tmp = zeros(Float32, (block_size,block_size,block_size))

    #Only calculate lower half of K3
    for i in 1:n_blocks_per_dim
        for j in 1:i
            for k in 1:j
                dim1_range = (block_size*(i-1) + 1):(block_size*i)
                dim2_range = (block_size*(j-1) + 1):(block_size*j)
                dim3_range = (block_size*(k-1) + 1):(block_size*k)
                K3_GPU_block =  mcc3_tensor(Ψ, phi, [dim1_range, dim2_range, dim3_range])
                copyto!(tmp, K3_GPU_block)

                #Cant copy off GPU directly into slice :/
                K3[dim1_range, dim2_range, dim3_range] = copyto!(K3[dim1_range, dim2_range, dim3_range],tmp)
            end
        end
    end

    for i in eachindex(K3)
        if abs(K3[i]) < tol
            K3[i] = 0.0
        end
    end

    return K3
end


function mcc3_tensor(Ψ::CuArray{Float32, 3}, phi::CuArray{Float32, 2}, idx_ranges::Vector{UnitRange{Int64}})

    block_sizes = length.(idx_ranges)
    K3 = CUDA.zeros(block_sizes...);

    phi_block1 = view(phi, :, idx_ranges[1])
    phi_block2 = view(phi, :, idx_ranges[2])
    phi_block3 = view(phi, :, idx_ranges[3])

    @tensor begin
        K3[n,m,l] = Ψ[i,j,k] * phi_block1[i,n] * phi_block2[j,m] * phi_block3[k,l]
    end

    return K3
end

"""
Converts third order forces constants, `Ψ` into third order modal coupling constants (MCC).
Does not divide MCC calculation into smaller chunks. This might exhaust GPU memory.
"""
function mcc3_tensor(Ψ::CuArray{Float32, 3}, phi::CuArray{Float32, 2})

    K3 = CUDA.zeros(size(Ψ));

    @tensor begin
        K3[n,m,l] = Ψ[i,j,k] * phi[i,n] * phi[j,m] * phi[k,l]
    end

    return K3
end
