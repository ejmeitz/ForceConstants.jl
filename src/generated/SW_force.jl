
# Energy functions used to generate force 
function Φ₂(A, B, ϵ, σ, p, q, a, rᵢⱼ)
    return A*ϵ*(B*((σ/rᵢⱼ)^p) - ((σ/rᵢⱼ)^q)) * exp(σ/(rᵢⱼ-(a*σ)))
end

function pair_potential_nounits(pot::StillingerWeberSilicon, dist_ij)   
    return Φ₂(pot.params.A, pot.params.B, ustrip(pot.params.ϵ), ustrip(pot.params.σ),
                 pot.params.p, pot.params.q, pot.params.a, dist_ij)
end

#Some simplifications are possible when using silicon params
function Φ₃_si(lambda_epsilon, gamma_sigma, cosθ₀, cosθᵢⱼₖ, r_cut, rᵢⱼ, rᵢₖ)
    return lambda_epsilon*((cosθᵢⱼₖ - cosθ₀)^2)*exp(gamma_sigma/(rᵢⱼ-r_cut))*exp(gamma_sigma/(rᵢₖ-r_cut))
end 

function three_body_potential_nounits(pot::StillingerWeberSilicon, rᵢⱼ, rᵢₖ, dist_ij, dist_ik)
    
    cosθᵢⱼₖ = dot(rᵢⱼ, rᵢₖ) / (dist_ij * dist_ik)

    return Φ₃_si(ustrip(pot.lambda_epsilon), ustrip(pot.gamma_sigma), pot.params.cosθ₀, cosθᵢⱼₖ, ustrip(pot.r_cut), dist_ij, dist_ik)
end

function generate_SW_force(path, pot::StillingerWeberSilicon)

    ri_vars = make_variables(:ri, 3)
    rj_vars = make_variables(:rj, 3)
    rk_vars = make_variables(:rk, 3)

    r_ij_vars = ri_vars .- rj_vars
    r_ik_vars = ri_vars .- rk_vars
    dist_ij = sqrt(sum(x -> x^2, r_ij_vars))
    dist_ik = sqrt(sum(x -> x^2, r_ik_vars))

    two_body_symb = pair_potential_nounits(pot, dist_ij)
    two_body_force = -FastDifferentiation.jacobian([two_body_symb], ri_vars)

    three_body_symb = three_body_potential_nounits(pot, r_ij_vars, r_ik_vars, dist_ij, dist_ik)
    three_body_force = -FastDifferentiation.jacobian([three_body_symb], ri_vars)

    force = two_body_force + three_body_force
    force_expr = make_Expr([force], [ri_vars; rj_vars; rk_vars], false, true)

    write(joinpath(path, "SW_force_test.jl"), string(force_expr))

end


#Symbolics.jl Version (won't take derivatives of array values)
# function generate_SW_force(path)

#     @variables ri[1:3] rj[1:3] rk[1:3]
#     @variables A B ϵ σ p q a
#     @variables lambda_epsilon gamma_sigma cosθ₀ cosθᵢⱼₖ

#     rij = ri .- rj
#     rik = ri .- rk
#     dist_ij = norm(rij)
#     dist_ik = norm(rik)
#     cosθᵢⱼₖ = dot(rij, rik) / (dist_ij * dist_ik)

#     two_body_expr = Φ₂(A, B, ϵ, σ, p, q, a, dist_ij)
#     three_body_expr = Φ₃(lambda_epsilon, gamma_sigma, cosθ₀, cosθᵢⱼₖ, σ*a, dist_ij, dist_ik)

#     D_ia = Differential(ri[1])
#     D_ib = Differential(ri[2])
#     D_ic = Differential(ri[3])

#     F_ia = expand_derivatives(-D_ia(two_body_expr + three_body_expr), true)
#     F_ib = expand_derivatives(-D_ib(two_body_expr + three_body_expr), true)
#     F_ic = expand_derivatives(-D_ic(two_body_expr + three_body_expr), true)
    
#     write(joinpath(path, "F_ia.jl"), string(build_function(F_ia, ri, rj, rk, cosθᵢⱼₖ, A, B, ϵ, σ, p, q, a, lambda_epsilon, gamma_sigma, cosθ₀)))
#     write(joinpath(path, "F_ib.jl"), string(build_function(F_ib, ri, rj, rk, cosθᵢⱼₖ, A, B, ϵ, σ, p, q, a, lambda_epsilon, gamma_sigma, cosθ₀)))
#     write(joinpath(path, "F_ic.jl"), string(build_function(F_ic, ri, rj, rk, cosθᵢⱼₖ, A, B, ϵ, σ, p, q, a, lambda_epsilon, gamma_sigma, cosθ₀)))
# end