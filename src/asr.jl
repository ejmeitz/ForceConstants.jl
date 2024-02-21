export ASR!

"""
Replaces entries in `ifc3` that could (and should) be calculated with the acoustic sum rule.
In ASR you fix (i,j) and (α, β, γ) and the sum over all k should be 0. Therefore, we 
can pick which k to solve for. This code chooses the ith element so that the self terms 
are calculated.
"""
function ASR!(ifc3::Array{T,3}, N_atoms, D) where T
    #Loop all self terms
    Threads.@threads for i in range(1,N_atoms)
        for j in range(1,N_atoms)
            #Fix row 
            for α in range(1,D)
                for β in range(1,D)
                    for γ in range(1,D)
                        ii_self = D*(i-1) + α; jj_self = D*(j-1) + β; kk_self = D*(i-1) + γ 

                        # Reset value and loop over row
                        ifc3[ii_self,jj_self,kk_self] = zero(T)
                        for k in range(1,N_atoms)
                            if k != i
                                kk = D*(k-1) + γ
                                ifc3[ii_self,jj_self,kk_self] -= ifc3[ii_self,jj_self,kk] 
                            end
                        end

                    end
                end
            end
        end
    end

    return ifc3
end

# function ASR!(ifc3::Array{T,3}, N_atoms, D) where T
#     #Loop all self terms
#     Threads.@threads for i in range(1,N_atoms)
#         for α in range(1,D)
#             for β in range(1,D)
#                 for γ in range(1,D)
#                     ii_self = D*(i-1) + α; jj_self = D*(i-1) + β; kk_self = D*(i-1) + γ 

#                     # Loop atoms k with this α, β, γ
#                     ifc3[ii_self,jj_self,kk_self] = zero(T)
#                     for k in range(1,N_atoms)
#                         if k != i
#                             ii = D*(i-1) + α; jj = D*(i-1) + β; kk = D*(k-1) + γ
#                             ifc3[ii_self,jj_self,kk_self] -= ifc3[ii,jj,kk] 
#                         end
#                     end

#                 end
#             end
#         end
#     end

#     return ifc3
# end