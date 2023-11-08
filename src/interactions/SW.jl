export StillingerWeberSilicon

struct SW_Params{E,S}
    ϵ::E #ep has energy units
    σ::S #simga has length units
    a::Float64 #cutoff = a*sigma
    λ::Float64
    γ::Float64
    cosθ₀::Float64
    A::Float64
    B::Float64
    p::Float64
    q::Float64
end

function SW_Params(ϵ, σ, a, λ, γ, cosθ₀, A, B , p, q)
    return SW_Params{typeof(ϵ), typeof(σ)}(ϵ, σ, a, λ, γ, cosθ₀,
                                             A, B , p, q)
end

struct StillingerWeberSilicon{R,L,E} <: ThreeBodyPotential
    r_cut::R
    length_unit::L
    energy_unit::E
    params::SW_Params
end


function StillingerWeberSilicon()
    #From LAMMPS Si.sw file
    sws_params = SW_Params(2.1683u"eV", 2.0951u"Å", 1.80, 21.0, 1.20,
         -1/3, 7.049556277,  0.6022245584,  4.0,  0.0)
    
    r_cut = sws_params.a*sws_params.σ
    length_unit = unit(sws_params.σ)
    energy_unit = unit(sws_params.ϵ)

    return StillingerWeberSilicon{typeof(r_cut),typeof(length_unit),typeof(energy_unit)}(
                r_cut, length_unit, energy_unit, sws_params)
end


function pair_potential(pot::StillingerWeberSilicon, rᵢⱼ::Vector)   
    return Φ₂(pot.params.A, pot.params.B, pot.params.ϵ, pot.params.σ, pot.params.p, pot.params.q, pot.params.a, norm(rᵢⱼ))
end

function pair_potential_nounits(pot::StillingerWeberSilicon, rᵢⱼ::Vector)   
    return Φ₂(pot.params.A, pot.params.B, ustrip(pot.params.ϵ), ustrip(pot.params.σ),
                 pot.params.p, pot.params.q, pot.params.a, norm(rᵢⱼ))
end

# function three_body_potential(pot::StillingerWeberSilicon, rᵢⱼ, rᵢₖ, rⱼᵢ, rⱼₖ, rₖᵢ, rₖⱼ)
#     dist_ij = norm(rᵢⱼ); dist_ik = norm(rᵢₖ)
#     dist_ji = norm(rⱼᵢ); dist_jk = norm(rⱼₖ)
#     dist_ki = norm(rₖᵢ); dist_kj = norm(rₖⱼ)
#     cosθᵢⱼₖ = dot(rᵢⱼ, rᵢₖ) / (dist_ij * dist_ik)
#     cosθⱼᵢₖ = dot(rⱼᵢ, rⱼₖ) / (dist_ji * dist_jk)
#     cosθₖᵢⱼ = dot(rₖᵢ, rₖⱼ) / (dist_ki * dist_kj)
#     h1 = Φ₃(pot.params.γ, pot.params.ϵ, cosθᵢⱼₖ, pot.params.cosθ₀, pot.params.γ, pot.params.γ,
#                  pot.params.σ, pot.params.σ, pot.params.a, pot.params.a, dist_ij, dist_ik)
#     h2 = Φ₃(pot.params.γ, pot.params.ϵ, cosθⱼᵢₖ, pot.params.cosθ₀, pot.params.γ, pot.params.γ,
#         pot.params.σ, pot.params.σ, pot.params.a, pot.params.a, dist_ji, dist_jk)
#     h3 = Φ₃(pot.params.γ, pot.params.ϵ, cosθₖᵢⱼ, pot.params.cosθ₀, pot.params.γ, pot.params.γ,
#         pot.params.σ, pot.params.σ, pot.params.a, pot.params.a, dist_ki, dist_kj)
#     # println("dist_ij $(dist_ij) dist_ik $(dist_ik) dist_ji $(dist_ji) dist_jk $(dist_jk) dist_ki $(dist_ki) dist_kj $(dist_kj)")
#     # println("h1: $(ustrip(h1)) h2 $(ustrip(h2)) h3 $(ustrip(h3))")
#     # println("cos1: $(ustrip(cosθᵢⱼₖ)) cos2 $(ustrip(cosθⱼᵢₖ)) cos3 $(ustrip(cosθₖᵢⱼ))")
#     return h1 + h2 + h3
# end
function three_body_potential(pot::StillingerWeberSilicon, rᵢⱼ, rᵢₖ)
    dist_ij = norm(rᵢⱼ); dist_ik = norm(rᵢₖ)
  
    cosθᵢⱼₖ = dot(rᵢⱼ, rᵢₖ) / (dist_ij * dist_ik)

    return Φ₃(pot.params.γ, pot.params.ϵ, cosθᵢⱼₖ, pot.params.cosθ₀, pot.params.γ, pot.params.γ,
                 pot.params.σ, pot.params.σ, pot.params.a, pot.params.a, dist_ij, dist_ik)

end

function three_body_potential_nounits(pot::StillingerWeberSilicon, rᵢⱼ, rᵢₖ, rⱼᵢ, rⱼₖ, rₖᵢ, rₖⱼ)
    dist_ij = norm(rᵢⱼ); dist_ik = norm(rᵢₖ)
    dist_ji = norm(rⱼᵢ); dist_jk = norm(rⱼₖ)
    dist_ki = norm(rₖᵢ); dist_kj = norm(rₖⱼ)
    cosθᵢⱼₖ = dot(rᵢⱼ, rᵢₖ) / (dist_ij * dist_ik)
    cosθⱼᵢₖ = dot(rⱼᵢ, rⱼₖ) / (dist_ji * dist_jk)
    cosθₖᵢⱼ = dot(rₖᵢ, rₖⱼ) / (dist_ki * dist_kj)
    h1 = Φ₃(pot.params.γ, ustrip(pot.params.ϵ), cosθᵢⱼₖ, pot.params.cosθ₀, pot.params.γ, pot.params.γ,
                 ustrip(pot.params.σ), ustrip(pot.params.σ), pot.params.a, pot.params.a, dist_ij, dist_ik)
    h2 = Φ₃(pot.params.γ, ustrip(pot.params.ϵ), cosθⱼᵢₖ, pot.params.cosθ₀, pot.params.γ, pot.params.γ,
        ustrip(pot.params.σ), ustrip(pot.params.σ), pot.params.a, pot.params.a, dist_ji, dist_jk)
    h3 = Φ₃(pot.params.γ, ustrip(pot.params.ϵ), cosθₖᵢⱼ, pot.params.cosθ₀, pot.params.γ, pot.params.γ,
        ustrip(pot.params.σ), ustrip(pot.params.σ), pot.params.a, pot.params.a, dist_ki, dist_kj)
    return h1 + h2 + h3
end
 
function Φ₂(A, B, ϵ, σ, p, q, a, rᵢⱼ)
    return A*ϵ*(B*((σ/rᵢⱼ)^p) - ((σ/rᵢⱼ)^q)) * exp(σ/(rᵢⱼ-(a*σ)))
end

function Φ₃(λᵢⱼₖ, ϵᵢⱼₖ, cosθᵢⱼₖ, cosθ₀, γᵢⱼ, γᵢₖ, σᵢⱼ, σᵢₖ, aᵢⱼ, aᵢₖ, rᵢⱼ, rᵢₖ)
    # println("term1: $(ustrip(λᵢⱼₖ*ϵᵢⱼₖ*((cosθᵢⱼₖ - cosθ₀)^2))) term2 $(ustrip(exp((γᵢⱼ*σᵢⱼ)/(rᵢⱼ -(aᵢⱼ*σᵢⱼ))))) term3 $(ustrip(exp((γᵢₖ*σᵢₖ)/(rᵢₖ-(aᵢₖ*σᵢₖ)))))")
    return λᵢⱼₖ*ϵᵢⱼₖ*((cosθᵢⱼₖ - cosθ₀)^2)*exp((γᵢⱼ*σᵢⱼ)/(rᵢⱼ -(aᵢⱼ*σᵢⱼ)))*exp((γᵢₖ*σᵢₖ)/(rᵢₖ-(aᵢₖ*σᵢₖ)))
end 
