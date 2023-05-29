
# function third_order(sys::SuperCellSystem{D}, pot::PairPotential) where D
#     N_atoms = n_atoms(sys)
#     Ψ = zeros(D*N_atoms,D*N_atoms,D*N_atoms)

#     unique_indices = with_replacement_combinations((1:N_atoms), 3)
#     for idx = unique_indicies
#         i,j,k = idx
#         for α in range(1,D)
#             for β in range(1,D)
#                 for γ in range(1,D)
#                     ii = 3*(i-1) + α
#                     jj = 3*(j-1) + β
#                     kk = 3*(k-1) + γ

#                     val = ϕ₃(pot, r_norm, rᵢⱼ, α, β, γ)
#                     Ψ[ii,jj,kk] = 1
#                 end
#             end
#         end
#     end
# end

# #Kronicker Delta
# δ(x,y) = ==(x,y)

# function ϕ₃(pot::PairPotential, r_norm, rᵢⱼ, α, β, γ)
#     Φ′ = potential_first_deriv(pot, r_norm)
#     Φ′′ = potential_second_deriv(pot, r_norm)
#     Φ′′′ = potential_third_deriv(pot, r_norm)

#     return (rᵢⱼ[α]*rᵢⱼ[β]*rᵢⱼ[γ]/(r_norm^3))*(Φ′′′ - (3*Φ′′/r_norm) + (3*Φ′/(r_norm^2))) +
#         ((rᵢⱼ[α]*δ(β,γ) + rᵢⱼ[β]*δ(γ,α) + rᵢⱼ[γ]*δ(α,β))/(r_norm^2))*(Φ′′ - (Φ′/r_norm))

# end