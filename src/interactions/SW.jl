#As one way, we assume that 
# lambda_ijk is equal to lambda_ik and eps_ijk is 
# equal to sqrt(lambda_ij*eps_ij*lambda_ik*eps_ik)/lambda_ik, 
# and all other parameters in the ijk line are for ik.

#https://www.afs.enea.it/software/lammps/doc17/html/pair_sw.html

struct SW_Params
    ϵ#ep has energy units
    σ #simga has length units
    a #cutoff = a*sigma
    λ
    γ
    cosθ₀
    A
    B
    p
    q
    r_cut
end

function SW_Params(ϵ, σ, a, λ, γ, cosθ₀, A, B , p, q)
    return SW_params(ϵ, σ, a, λ, γ, cosθ₀, A, B , p, q, a*σ)
end

struct StillingerWeberSilicon{L,E,R} <: ThreeBodyPotential
    length_unit::L
    energy_unit::E
    params::SW_Params #Dict{NTuple{3,Symbol},SW_Params}
end

function StillingerWeberSilicon(params::SW_Params)
    @assert unit(params.σ_ij) == unit(params.σ_ik) "Mismatched length units in SW"
    @assert unit(params.ϵ_ij) == unit(params.ϵ_ijk) "Mismatched energy units in SW"
    # @assert length(params ∈ [1,8,27]) "Invalid number of parameters"
    # @assert validate_permutations(params) "Missing element permutation"

    #Pull out units, assume they use consistent units
    len_unit = unit(first(params).second.σ_ij)
    eng_unit = unit(first(params).second.ϵ_ij)

    #Pull out values where last 2 are same (used for pair evals per LAMMPS docs)
    # two_body_inters = {k => v for (k,v) in params if all(k[2] == k[3])}

    return SW{typeof(len_unit),typeof(eng_unit),length(params)}(len_unit, eng_unit, params)
end

function validate_permutations(params)
    if length(params) == 1
        elements = keys(params[1])
        return all(elements[1] .== elements) #all must be the same
    elseif length(params) == 8
        # Check only 2 unique elements
        # Check there is all possible combos
        return only_2_elements && all_perms_present
    elseif length(params) == 27
        #Check for 3 unique elements
        #Check there is all possible combos
        return only_3_elements && all_perms_present
    else
        throw(error("Invalid number of params"))
    end
end


function potential(pot::StillingerWeberSilicon, rᵢ, rⱼ, rₖ, e1, e2, e3)
    params = pot.params[(e1,e2,e3)]
    rᵢⱼ = rᵢ .- rⱼ
    rᵢₖ = rᵢ .- rₖ
    θᵢⱼₖ = dot(rᵢⱼ, rᵢₖ) / (norm(rᵢⱼ) * norm(rᵢₖ))

    λᵢⱼₖ = pot.ϵ
    # sqrt(lambda_ij*eps_ij*lambda_ik*eps_ik)/lambda_ik
    return Φ₂(params.A, params.B, params.ϵ, params.p, params.q, params.a, rᵢⱼ) +
                  Φ₃(, rᵢⱼ, rᵢₖ, θᵢⱼₖ)

end

function Φ₂(A, B, ϵ, p, q, a, rᵢⱼ)
    return A*ϵ*(B*((σ/rᵢⱼ)^p) - ((σ/rᵢⱼ)^q)) * exp(σ/(rᵢⱼ-(a*σ)))
end

function Φ₃(λᵢⱼₖ, ϵᵢⱼₖ, cosθ₀, γᵢⱼ, γᵢₖ, σᵢⱼ, σᵢₖ, aᵢⱼ, aᵢₖ, rᵢⱼ, rᵢₖ, θᵢⱼₖ)
    return λᵢⱼₖ*ϵᵢⱼₖ*((cos(θᵢⱼₖ) - cosθ₀)^2)*exp((γᵢⱼ*σᵢⱼ)/(rᵢⱼ -(aᵢⱼ*σᵢⱼ)))*exp((γᵢₖ*σᵢₖ)/(rᵢₖ-(aᵢₖ*σᵢₖ)))
end

function pair_force()

end

function three_body_force()

end