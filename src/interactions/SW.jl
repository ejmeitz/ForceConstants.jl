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

struct StillingerWeberSilicon{R,GS,LE} <: ThreeBodyPotential
    r_cut::R
    params::SW_Params
    gamma_sigma::GS
    lambda_epsilon::LE
end


function StillingerWeberSilicon()
    #From LAMMPS Si.sw file
    sws_params = SW_Params(2.1683u"eV", 2.0951u"Å", 1.80, 21.0, 1.20,
        -0.333333333333, 7.049556277,  0.6022245584,  4.0,  0.0)
    
    r_cut = sws_params.a*sws_params.σ
    gamma_sigma = sws_params.σ*sws_params.γ
    lambda_epsilon = sws_params.ϵ*sws_params.λ

    return StillingerWeberSilicon{typeof(r_cut), typeof(gamma_sigma), typeof(lambda_epsilon)}(
                r_cut, sws_params, gamma_sigma, lambda_epsilon)
end

energy_unit(pot::StillingerWeberSilicon) = unit(pot.params.ϵ)
length_unit(pot::StillingerWeberSilicon) = unit(pot.params.σ)


function pair_potential(pot::StillingerWeberSilicon, dist_ij)   
    return Φ₂(pot.params.A, pot.params.B, pot.params.ϵ, pot.params.σ, pot.params.p, pot.params.q, pot.params.a, dist_ij)
end

function pair_potential_nounits(pot::StillingerWeberSilicon, dist_ij)   
    return Φ₂(pot.params.A, pot.params.B, ustrip(pot.params.ϵ), ustrip(pot.params.σ),
                 pot.params.p, pot.params.q, pot.params.a, dist_ij)
end


function three_body_potential(pot::StillingerWeberSilicon, rᵢⱼ, rᵢₖ, dist_ij, dist_ik)
    
    cosθᵢⱼₖ = dot(rᵢⱼ, rᵢₖ) / (dist_ij * dist_ik)

    return Φ₃_si(pot.lambda_epsilon, pot. gamma_sigma, pot.params.cosθ₀, cosθᵢⱼₖ, pot.r_cut, dist_ij, dist_ik)

end

function three_body_potential_nounits(pot::StillingerWeberSilicon, rᵢⱼ, rᵢₖ, dist_ij, dist_ik)
    
    cosθᵢⱼₖ = dot(rᵢⱼ, rᵢₖ) / (dist_ij * dist_ik)

    return Φ₃_si(ustrip(pot.lambda_epsilon), ustrip(pot.gamma_sigma), pot.params.cosθ₀, cosθᵢⱼₖ, ustrip(pot.r_cut), dist_ij, dist_ik)
end


 
function Φ₂(A, B, ϵ, σ, p, q, a, rᵢⱼ)
    return A*ϵ*(B*((σ/rᵢⱼ)^p) - ((σ/rᵢⱼ)^q)) * exp(σ/(rᵢⱼ-(a*σ)))
end

#Some simplifications are possible when using silicon params
function Φ₃_si(lambda_epsilon, gamma_sigma, cosθ₀, cosθᵢⱼₖ, r_cut, rᵢⱼ, rᵢₖ)
    return lambda_epsilon*((cosθᵢⱼₖ - cosθ₀)^2)*exp(gamma_sigma/(rᵢⱼ-r_cut))*exp(gamma_sigma/(rᵢₖ-r_cut))
end 

# function Φ₃(λᵢⱼₖ, ϵᵢⱼₖ, cosθᵢⱼₖ, cosθ₀, γᵢⱼ, γᵢₖ, σᵢⱼ, σᵢₖ, aᵢⱼ, aᵢₖ, rᵢⱼ, rᵢₖ)
#     return λᵢⱼₖ*ϵᵢⱼₖ*((cosθᵢⱼₖ - cosθ₀)^2)*exp((γᵢⱼ*σᵢⱼ)/(rᵢⱼ-(aᵢⱼ*σᵢⱼ)))*exp((γᵢₖ*σᵢₖ)/(rᵢₖ-(aᵢₖ*σᵢₖ)))
# end 

