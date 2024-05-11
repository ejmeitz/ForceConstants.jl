


function get_ALAMODE_data(ifc2::Matrix, primitive_atom_idxs, in_units::Symbol; D = 3)

    alm_ifc2 = Dict{Tuple{Int, Int}, eltype(ifc2)}()
    # cart_map = Dict(1 => "x", 2 => "y", 3 => "z")
    N_atoms = Int64(size(ifc2)[1] / D)

    if in_units == :METAL
        conv = (1.889725988*1.889725988)*13.60569301
    elseif in_units == :REAL
        conv = (1.889725988*1.889725988)*313.75470835207074
    else
        error("Unknown unit system, must be METAL or REAL")
    end

    for i in primitive_atom_idxs
        for j in 1:N_atoms
            for α in 1:D
                for β in 1:D
                    # key = ("$(i)$(cart_map[α])", "$(j)$(cart_map[β])")
                    ii = D*(i-1)+α; jj = D*(j-1)+β
                    if ifc2[ii,jj] != 0.0
                        alm_ifc2[(ii-1,jj-1)] = ifc2[ii,jj]/conv
                    end
                end
            end
        end
    end

    return alm_ifc2
end

function get_ALAMODE_data(ifc3::Array{T,3}, primitive_atom_idxs, in_units; D = 3) where T

    alm_ifc3 = Dict{Tuple{Int, Int, Int}, eltype(ifc3)}()
    N_atoms = Int64(size(ifc3)[1] / D)

    if in_units == :METAL
        conv = (1.889725988*1.889725988*1.889725988)*13.60569301
    elseif in_units == :REAL
        conv = (1.889725988*1.889725988*1.889725988)*313.75470835207074
    else
        error("Unknown unit system, must be METAL or REAL")
    end

    for i in primitive_atom_idxs
        for j in 1:N_atoms
            for k in 1:N_atoms
                for α in 1:D
                    for β in 1:D
                        for γ in 1:D
                            ii = D*(i-1)+α; jj = D*(j-1)+β; kk = D*(k-1)+γ
                            if ifc3[ii,jj,kk] != 0.0
                                alm_ifc3[(ii-1,jj-1,kk-1)] = ifc3[ii,jj,kk]/conv
                            end
                        end
                    end
                end
            end
        end
    end

    return alm_ifc3
end


function write_alm_pairs(ifc_kv::AbstractDict, outpath)
    open(outpath, "w") do f
        k = keys(ifc_kv)
        D = length(first(k))
        idxs = [[key[i] for key in k] for i in 1:D]
        vals = values(ifc_kv)
        k_out = zeros(eltype(first(k)), length(k), D)
        v_out = collect(vals)
        for i in 1:length(k)
            for j in 1:D
                k_out[i,j] = idxs[j][i]
            end
        end
        writedlm(f, Any[k_out v_out])
    end
end

# Averages over non_zero_idxs indexes
# Sets zero_idxs to zero
#Applies ASR at end
function clean_avg_ifc!(avg_ifc::AbstractArray, ifc_0K::AbstractArray; digits = 8)

    N_modes = size(avg_ifc)[1]
    N_atoms = Int64(N_modes / 3)

    ifc_0K = round.(ifc_0K, digits = digits)
    ifc_0K_unique = unique(ifc_0K)
    filter!( x -> x != 0.0, ifc_0K_unique)

    for i in eachindex(ifc_0K_unique)
        non_zero_idx = findall(x -> x == ifc_0K_unique[i], ifc_0K)
        m = mean(avg_ifc[non_zero_idx])
        avg_ifc[non_zero_idx] .= m
    end

    zero_idxs = findall(x -> x == 0.0, ifc_0K)
    avg_ifc[zero_idxs] .= 0.0

    ASR!(avg_ifc, N_atoms, 3)  

    return avg_ifc
end

function ASR_test!(ifc2::Matrix{T}, N_atoms, D) where T

    Threads.@threads for i in range(1, N_atoms) # index of block matrix
        for α in range(1,D)
            for β in range(α,D)
                ii = D*(i-1) + α
                jj = D*(i-1) + β # i == j because we're on diagonal
                ifc2[ii,jj] = zero(T)
                ifc2[ii,jj] = -1*sum(ifc2[ii, β:D:end])
                ifc2[jj,ii] = ifc2[ii,jj]
            end
        end
    end

    return ifc2
end

function ASR_test!(ifc3::Array{T,3}, N_atoms, D) where T
    #Loop all self terms
    Threads.@threads for i in range(1,N_atoms)
        for j in range(1,N_atoms)
            for α in range(1,D)
                ii_self = D*(i-1) + α
                for β in range(α,D)
                    jj_self = D*(j-1) + β
                    for γ in range(β,D)
                        kk_self = D*(i-1) + γ 

                        ifc3[ii_self,jj_self,kk_self] = zero(T)
                        ifc3[ii_self,jj_self,kk_self] = -sum(ifc3[ii_self,jj_self,γ:D:end]) 

                        ifc3[ii_self,kk_self,jj_self] = ifc3[ii_self,jj_self,kk_self]
                        ifc3[jj_self,ii_self,kk_self] = ifc3[ii_self,jj_self,kk_self]
                        ifc3[jj_self,kk_self,ii_self] = ifc3[ii_self,jj_self,kk_self]
                        ifc3[kk_self,ii_self,jj_self] = ifc3[ii_self,jj_self,kk_self]
                        ifc3[kk_self,jj_self,ii_self] = ifc3[ii_self,jj_self,kk_self]
                    end
                end
            end
        end
    end

    return ifc3
end


function clean_avg_ifc_safe!(avg_ifc::AbstractMatrix, sys::SuperCellSystem, r_cut, tol)
    N_atoms = n_atoms(sys)

    rᵢⱼ = similar(sys.atoms.position[1])
    r_cut_sq = r_cut*r_cut
    for i in 1:N_atoms
        for j in i+1:N_atoms
            rᵢⱼ .= sys.atoms.position[i] .- sys.atoms.position[j]
            nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
            dist_ij_sq = sum(x -> x^2, rᵢⱼ)

            if dist_ij_sq > r_cut_sq
                avg_ifc[3*(i-1)+1 : 3*i, 3*(j-1)+1 : 3*j] .= 0.0
                avg_ifc[3*(j-1)+1 : 3*j, 3*(i-1)+1 : 3*i] .= 0.0
            end
        end
    end

    avg_ifc[abs.(avg_ifc) .< tol] .= 0.0

    ASR_test!(avg_ifc, N_atoms, 3)
    

    return avg_ifc

end

function clean_avg_ifc_safe!(avg_ifc::Array{T,3}, sys::SuperCellSystem, r_cut, tol) where T
    N_atoms = n_atoms(sys)

    rᵢⱼ = similar(sys.atoms.position[1])
    r_cut_sq = r_cut*r_cut
    for i in 1:N_atoms
        for j in 1:N_atoms
            if i != j
                rᵢⱼ .= sys.atoms.position[i] .- sys.atoms.position[j]
                nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
                dist_ij_sq = sum(x -> x^2, rᵢⱼ)

                if dist_ij_sq < r_cut_sq
                    for k in range(j+1, N_atoms)
                        if i != k
                            rᵢₖ .= sys.atoms.position[i] .- sys.atoms.position[k]
                            nearest_mirror!(rᵢₖ, sys.box_sizes_SC)
                            dist_ik_sq = sum(x -> x^2, rᵢₖ)
                            if dist_ik_sq > r_cut_sq
                                # SET ijk BLOCK TO ZERO AND ALL SYMMETRIES
                            end
                        end
                    end
                else
                    # SET WHOLE k SLICE TO 0 AND ALL SYMMETRIES
                end


            end
        end
    end

    avg_ifc[abs.(avg_ifc) .< tol] .= 0.0

    ASR_test!(avg_ifc, N_atoms, 3)
    

    return avg_ifc

end


function clean_avg_ifc_by_uc(sys, eq_posns, primitive_atoms, equivalent_atoms, r_cut, avg_ifc)

    # Find interactions possible for atoms in conventional cell
    r_cut_sq = r_cut*r_cut
    interactions = Dict(pa => [] for pa in primitive_atoms)
    r_ij = similar(eq_posns[1])
    for i in primitive_atoms
        pa_atom = eq_posns[i]
        for (j,eq) in enumerate(eq_posns)
            if i != j
                r_ij .= eq .- pa_atom
                nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
                dist_ij_sq = sum(x -> x^2, r_ij)
                if dist_ij_sq < r_cut_sq
                    push!(interactions[pa], eq)
                end
            end
        end
    end

    # Loop all interactions in super cell and figure out which
    # primitive cell interaction they map to
    for i in 1:n_atoms(sys)
        # Map this atom back to primitive_cell
        cc_atom_i = equivalent_atoms[i] 

        if cell_idx != [0,0,0] 
            for j in (i+1):n_atoms(sys)

                # Map one atom back to its equivalent atom in unit cell
            

            end
        end
    # Average interactions in supercell that are the same


end

# out_path(T) = "Z:/emeitz/Data/ForceConstants/AvgINM_SW/AvgIFC_SW_3UC_$(T)K_CLEAN.jld2"
# for T in temps
#     @info "$T"
#     dynmat_T, ifc3_T = load(tep_path(T), "dynmat", "F3")
#     dynmat_T = clean_avg_ifc!(dynmat_T, dynmat_0K)
#     freqs_sq, phi = get_modes(dynmat_T)
#     ifc3_T = clean_avg_ifc!(ifc3_T, ifc3_0K)
#     K3 = Array(mcc3(CuArray{Float32}(ifc3_T ./ sqrt(m^3)), CuArray{Float32}(phi)))
#     jldsave(out_path(T), dynmat = dynmat_T, freqs_sq = freqs_sq, F3 = ifc3_T, K3 = K3, phi = phi)
# end
