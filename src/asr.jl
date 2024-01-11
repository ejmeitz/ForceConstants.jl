function ASR!(ifc3::Array{T,3}, N_atoms, D) where T
    #Loop all self terms
    Threads.@threads for i in range(1,N_atoms)
        for α in range(1,D)
            for β in range(1,D)
                for γ in range(1,D)
                    ii_self = D*(i-1) + α; jj_self = D*(i-1) + β; kk_self = D*(i-1) + γ 

                    # Loop atoms k with this α, β, γ
                    for k in range(1,N_atoms)
                        if k != i
                            ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(k-1) + γ
                            ifc3[ii_self,jj_self,kk_self] -= ifc3[ii,jj,kk] 
                        end
                    end

                end
            end
        end
    end

    return ifc3
end