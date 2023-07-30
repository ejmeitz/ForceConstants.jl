export mcc3


"""
Converts third order forces constants, `Ψ` into third order modal coupling constants (MCC).
Does not divide MCC calculation into smaller chunks. This might exhaust GPU memory.
"""
function mcc3(Ψ::CuArray{Float32, 3}, phi::CuArray{Float32, 2}, tol = 1f-12)

    K3 = CUDA.zeros(Float32, size(Ψ));

    @tensor begin
        K3[n,m,l] = Ψ[i,j,k] * phi[i,n] * phi[j,m] * phi[k,l]
    end

    K3_CPU = Array(K3)

    for i in eachindex(K3_CPU)
        if abs(K3_CPU[i]) < tol
            K3_CPU[i] = 0.0f32
        end
    end

    return K3_CPU
end

"""
Converts third order forces constants, `Ψ` into third order modal coupling constants (MCC). The
parameter `block_size` specifies problem size when calculating the MCC. For example, if 
`block_size` is 100, the block MCC will be calculated in 100x100x100 blocks to save GPU memory.
The `phi` matrix will be automatically truncated to adjust for this. Try to minimize block_size.
"""
function mcc3(Ψ::CuArray{Float32, 3}, phi::CuArray{Float32, 2}, block_size::Int; tol = 1.0f-12)

    @assert size(phi)[1] % block_size == 0
    @assert block_size > 0

    n_blocks_per_dim = Int(size(phi)[1] / block_size)

    #Keep large storage on CPU
    K3_CPU = zeros(Float32, size(Ψ))

    #Only calculate lower half of K3
    for i in 1:n_blocks_per_dim
        for j in 1:i
            for k in 1:j
                dim1_range = (block_size*(i-1) + 1):(block_size*i)
                dim2_range = (block_size*(j-1) + 1):(block_size*j)
                dim3_range = (block_size*(k-1) + 1):(block_size*k)
                K3_GPU_block =  mcc3_blocked(Ψ, phi, [dim1_range, dim2_range, dim3_range])
                copyto!(view(K3_CPU, dim1_range, dim2_range, dim3_range), K3_GPU_block)
            end
        end
    end

    #Enforce symmetry of K3
    throw(error("Not yet implemented"))
    for i in 1:n_blocks_per_dim
        for j in 1:i
            for k in 1:j
                idx_combos = multiset_permutations([i,j,k], 3)
                filter!(p -> p != [i,j,k] , collect(idx_combos))
                dim1_range_orig = (block_size*(i-1) + 1):(block_size*i)
                dim2_range_orig = (block_size*(j-1) + 1):(block_size*j)
                dim3_range_orig = (block_size*(k-1) + 1):(block_size*k)
                for idx in idx_combos
                    dim1_range = (block_size*(idx[1]-1) + 1):(block_size*idx[1])
                    dim2_range = (block_size*(idx[2]-1) + 1):(block_size*idx[2])
                    dim3_range = (block_size*(idx[3]-1) + 1):(block_size*idx[3])
                    copyto!(view(K3_CPU, dim1_range, dim2_range, dim3_range),
                            view(K3_CPU, dim1_range_orig, dim2_range_orig, dim3_range_orig))
                end
            end
        end
    end

    for i in eachindex(K3_CPU)
        if abs(K3_CPU[i]) < tol
            K3_CPU[i] = 0.0f32
        end
    end

    return K3_CPU
end


function mcc3_blocked(Ψ::CuArray{Float32, 3}, phi::CuArray{Float32, 2}, idx_ranges::Vector{UnitRange{Int64}})

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

