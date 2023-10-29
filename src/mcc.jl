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

function mcc4_stupid!(ifc4_mw::SparseForceConstants{4}, phi, tol; nthreads = Threads.nthreads())
   
    mcc_idxs = with_replacement_combinations(1:N_modes, 4) #a casual 14 billion for 4th order

    chunks = Iterators.partition(mcc_idxs, cld(length(mcc_idxs), nthreads))

    @assert length(chunks) == nthreads

    tasks = map(chunks) do chunk

        Base.Threads.@spawn begin
            χ_thd = FC_val{Float64, 4}[]
            
            #& do on GPU with mapreduce??
            for mcc_idx in chunk #* this is incredibly slow, indexes being Int16 helps
                mcc_val = 0.0
                phi_n = view(phi,:,mcc_idx[1]); phi_m = view(phi,:,mcc_idx[2])
                phi_l = view(phi,:,mcc_idx[3]); phi_o = view(phi,:,mcc_idx[4])
                @turbo for idx in eachindex(ifc4.values)
                    mcc_val += vals[idx] * phi_n[is[idx]] * phi_m[js[idx]] * phi_l[ks[idx]] * phi_o[os[idx]]
                end

                if abs(mcc_val) > tol
                    push!(χ_thd, mcc_val)
                end
            end

            χ_thd
        end
    end

    χ_sparse_arrays = fetch.(tasks)
    return reduce(vcat, χ_sparse_arrays) #combine into one array

end

#& re-write to take CuArray since we dont want to move to GPU every time we calculate
# function gpu_k4_kernel_SA(cuF4mw::StructArray{F4_val}, cuPhi1, cuPhi2, cuPhi3, cuPhi4)
#     #Define function for map-reduce
#     f = (f4_data) -> f4_data.v * cuPhi1[f4_data.i] * cuPhi2[f4_data.j] * cuPhi3[f4_data.k] * cuPhi4[f4_data.l]
#    return mapreduce(f, +, cuF4mw)
# end

# function gpu_k4_kernel(cuF4mw::CuArray{F4_val}, cuPhi1, cuPhi2, cuPhi3, cuPhi4)
#     #Define function for map-reduce
#     f = (f4_data) -> f4_data.v * cuPhi1[f4_data.i] * cuPhi2[f4_data.j] * cuPhi3[f4_data.k] * cuPhi4[f4_data.l]
#    return mapreduce(f, +, cuF4mw)
# end

# struct F4_val
#     v::Float32
#     i::Int32
#     j::Int32
#     k::Int32
#     l::Int32
# end


# N_modes = 768
# N = 20000000
# q = rand(Float32, N_modes)
# random_numbers = rand(Float32,N);
# random_idxs = Int32.(rand(1:N_modes,(N,4)))
# sparse_test = [F4_val(random_numbers[i], random_idxs[i,:]...) for i in 1:N]
# sparse_test_SA = StructArray{F4_val}(v = random_numbers, i = random_idxs[:,1], j = random_idxs[:,2], k = random_idxs[:,3], l = random_idxs[:,4])
# cuF4mw = replace_storage(CuArray, sparse_test_SA)
# phi = rand(N_modes, N_modes)
# cuPhi = CuArray{Float32}(phi)
# gpu_k4_kernel_SA(cuF4mw, view(cuPhi,:,13), view(cuPhi,:,131), view(cuPhi,:,245), view(cuPhi,:,613))
# cuF4_vec = CuArray(sparse_test)
# gpu_k4_kernel(cuF4_vec,view(cuPhi,:,13), view(cuPhi,:,131), view(cuPhi,:,245), view(cuPhi,:,613))

#################################

# using CUDA
# import CUDA: i32
# N = 10000000
# N_modes = 768;
# rand_nums = CUDA.rand(N);
# rand_idxs = rand(1:N_modes, (N,3));
# is = CuArray{Int32}(rand_idxs[:,1])
# js = CuArray{Int32}(rand_idxs[:,2])
# ks = CuArray{Int32}(rand_idxs[:,3])

# phi = CUDA.rand(N_modes, N_modes)
# mcc_all = zeros(Float32, (N_modes, N_modes, N_modes))

# #256 x 256 x N_modes
# mcc_chunk_size = 256
# mcc_chunk = CUDA.zeros(Float32, (mcc_chunk_size, mcc_chunk_size, N_modes))


# #CUDA THREAD DIMENSIONS
# const BLOCK_SIDE::Int32 = 32
# const BLOCK_SIZE::Int32 = BLOCK_SIDE*BLOCK_SIDE

# const GRID_SIDE::Int64 = Int64(ceil(mcc_chunk_size/BLOCK_SIDE))
# @cuda blocks = (GRID_SIDE,GRID_SIDE) threads = (BLOCK_SIDE,BLOCK_SIDE) mcc3_custom_kernel(Int32(N), Int32(0), Int32(N_modes), mcc_chunk, phi,rand_nums,is,js,ks)

# #Each thread calculates fiber of K_nml
# function mcc3_custom_kernel(num_psi_idxs::Int32, mcc_offset::Int32, N_modes::Int32, mcc_block::CuDeviceArray{Float32,3}, phi::CuDeviceArray{Float32,2},
#     Ψ_vals::CuDeviceArray{Float32,1}, is::CuDeviceArray{Int32,1}, js::CuDeviceArray{Int32,1}, ks::CuDeviceArray{Int32,1}) 

#     n::Int32 = (blockIdx().x-1i32) * blockDim().x + threadIdx().x
#     m::Int32 = (blockIdx().y-1i32) * blockDim().y + threadIdx().y
#     @cuprintln "$n $m"
#     linearThreadIdx::Int32 = (blockDim().x * (threadIdx().y-1i32)) + threadIdx().x #& OFF BY 1 iSSUES?

#     #To store rows of phi matrix needed by current values of Psi
#     phi_i = CuStaticSharedArray(Float32, BLOCK_SIZE)
#     phi_j = CuStaticSharedArray(Float32, BLOCK_SIZE)
#     phi_k = CuStaticSharedArray(Float32, BLOCK_SIZE)

#     # Treat thread index as index into Ψ now and accumulate values into MCC
#     for i in 1:BLOCK_SIZE:num_psi_idxs
#         psi_idx::Int32 = i + linearThreadIdx

#         if psi_idx <= num_psi_idxs
#             #Move Phi Data Into Shared Memory
#             phi_i[linearThreadIdx] = phi[is[psi_idx],linearThreadIdx] #* this read is not in order of memory
#             phi_j[linearThreadIdx] = phi[js[psi_idx],linearThreadIdx]
#             phi_k[linearThreadIdx] = phi[ks[psi_idx],linearThreadIdx] #& This might be the only one that needs to be in shared mem
#             CUDA.sync_threads()

#             #Calculate fiber of MCC
#             # if m >= n #* this could make it slower??
#             for l in 1:N_modes
#                 mcc_block[n,m,l] += Ψ_vals[psi_idx]*phi_i[n + mcc_offset]*phi_j[m + mcc_offset]*phi_k[l]
#             end
#             # end
#             CUDA.sync_threads()
#         end
#     end
#     return
# end
