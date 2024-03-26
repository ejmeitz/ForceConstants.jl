export ASR!, asr_satisfied


function ASR!(ifc3::Array{T,3}, N_atoms, D) where T
    #Loop all self terms
    Threads.@threads for i in range(1,N_atoms)
        for j in range(1,N_atoms)
            for α in range(1,D)
                ii_self = D*(i-1) + α
                for β in range(1,D)
                    jj_self = D*(j-1) + β
                    for γ in range(1,D)
                        kk_self = D*(i-1) + γ 

                        ifc3[ii_self,jj_self,kk_self] = zero(T)
                        ifc3[ii_self,jj_self,kk_self] = -sum(ifc3[ii_self,jj_self,γ:D:end]) 
                    end
                end
            end
        end
    end

    return ifc3
end


function asr_satisfied(ifc2::AbstractArray{T,2}, N_atoms, D, tol; verbose = false) where T

    satisfied = true
    Threads.@threads for i in range(1,N_atoms)
        for α in range(1,D)
            for β in range(1,D)
                s = sum(ifc2[D*(i-1) + α, β:D:end])
                if sum(ifc2[D*(i-1) + α, β:D:end]) > tol
                    satisfied = false
                    verbose && println("Set i = $i, α = $α, β = $β does not satisfy ASR. Sum was: $s")
                end
            end
        end
    end

    return satisfied

end


function asr_satisfied(ifc3::Array{T,3}, N_atoms, D, tol; verbose = false) where T
    satisfied = true
    for k in range(1,N_atoms)
        for γ in range(1,D)
            ifc3_view = view(ifc3, :, :, D*(k-1) + γ)
            verbose && println("=======\n On k = $k, γ = $γ \n=======")
            slice_satisfied = asr_satisfied(ifc3_view, N_atoms, D, tol, verbose = verbose)
            satisfied &= slice_satisfied
        end
    end
    return satisfied
end

function rotational_invariance_satisfied()

end

function check_huang()

end