# export TersoffSilicon

# struct TersoffParams{E,IL,L}
#     m
#     γ
#     λ₃::IL
#     c
#     d
#     cosθ₀
#     n
#     β
#     λ₂::IL
#     B::E
#     R::L
#     D::L
#     λ₁::IL
#     A::E
#     δ::L
# end

# function TersoffParams(m, γ, λ₃, c, d, cosθ₀, n, β, λ₂, B, R, D, λ₁, A)
#     return TersoffParams{typeof(B), typeof(λ₁), typeof(R)}(m, γ, λ₃, c, d, cosθ₀, n, β, λ₂, B, R, D, λ₁, A)
# end

# struct TersoffSilicon{R} <: ThreeBodyPotential
#     r_cut::R
#     params::TersoffParams
# end


# function TersoffSilicon(δ = 0.0u"Å")
#     #From LAMMPS Si.sw file
#     sws_params = SW_Params(3.0, 1.0, 1.3258u"Å^-1", 4.8381, 2.0417, 0.0,
#                  22.956, 0.33675, 1.3258u"Å^-1", 95.373u"eV", 
#                  3.0u"Å", 0.2u"Å", 3.2394u"Å^-1", 3264.7u"eV")
    
#     r_cut = sws_params.R + sys_params.D

#     return TersoffSilicon{typeof(r_cut)}(r_cut, sws_params, δ)
# end

# energy_unit(pot::TersoffSilicon) = unit(pot.params.A)
# length_unit(pot::TersoffSilicon) = unit(pot.r_cut)

# function g(θ, ts::TersoffSilicon) 
#     return ts.γ*(1 + (ts.c/ts.d)^2 - ((ts.c^2)/(ts.d^2 + (cos(θ) - ts.cosθ₀)^2)))
# end

# function fₐ(r, ts::TersoffSilicon)
#     return -ts.B.*exp.(-ts.λ₂*r)
# end

# function fᵣ(r, ts::TersoffSilicon)
#     return ts.A.*exp.(-ts.λ₁*r)
# end

# function bᵢⱼ(ζᵢⱼ, ts::TersoffSilicon)
#     return (1 + ((ts.β*ζᵢⱼ).^ts.n))^(-1/(2*ts.n))
# end

# function fc(r, ts::TersoffSilicon)
#     if r < ts.R - ts.D
#         return 1
#     elseif ts.R - ts.D < r && r < ts.R + ts.D
#         return 0.5*(1 - sin.(0.5*pi*(r - ts.R)/ts.D))
#     else
#         return 0
#     end
# end

# function Vᵢⱼ(rᵢⱼ, bᵢⱼ)
#     return fc(rᵢⱼ .+ ts.δ, ts)*(fᵣ(rᵢⱼ.+ ts.δ, ts) + bᵢⱼ*fₐ(rᵢⱼ.+ ts.δ,ts))
# end

# #& This is precomputed each iteration??
# function ζᵢⱼ(ts::TersoffSilicon, N_atoms::Int, r::AbstractVector{<:AbstractVector}, box_sizes, i::Int)
#     res = 0.0
#     r_cut_sq = ts.r_cut*ts.r_cut
#     rᵢⱼ = similar(r[1])
#     rᵢₖ = similar(r[1])

#     for j in range(1, N_atoms)
#         rᵢⱼ .= r[i] .- r[j]
#         nearest_mirror!(rᵢⱼ, box_sizes)
#         dist_ij_sq = dot(rᵢⱼ, rᵢⱼ)

#         if dist_ij_sq < r_cut_sq
#             for k in range(1,N_atoms)
#                 rᵢₖ .= r[i] .- r[k]
#                 nearest_mirror!(rᵢₖ, box_sizes)
#                 dist_ik_sq = dot(rᵢₖ, rᵢₖ)
#                 if i != k && dist_ik_sq < r_cut_sq
#                     θᵢⱼₖ = dot(rᵢⱼ, rᵢₖ) / sqrt(dist_ij_sq * dist_ik_sq)
#                     res += fc(rᵢₖ .+ δ, ts).*g(θᵢⱼₖ, ts).*exp.(λ₃^ts.m*(rᵢⱼ - rᵢₖ))^m
#                 end
#             end
#         end
#     end
# end


# function tersoff_energy_loop(pot::TersoffSilicon, sys::SuperCellSystem)

#     U_total = 0.0*energy_unit(pot)

#     # Calculate ζᵢⱼ and bᵢⱼ
#     #TODO

#     for i in range(1,N_atoms)
#         for j in range(i+1, N_atoms)             
#             r = posns[i] .- posns[j]
#             nearest_mirror!(r, box_sizes)
#             dist = norm(r)
#             if dist < r_cut
#                 U_total += potential(pot,dist)
#             end
#         end
#     end

#     return U_total



# end