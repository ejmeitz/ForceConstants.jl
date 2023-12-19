function get_three_body_derivs(pot::StillingerWeberSilicon)

    #2nd order part
    r_ij_vars = make_variables(:rij, D) #rᵢ - rⱼ
    r_ij_norm = sqrt(sum(x -> x^2, r_ij_vars))

    pot2_symbolic = pair_potential_nounits(pot, r_ij_norm)
    H2_symbolic = hessian(pot2_symbolic, r_ij_vars)
    H2_exec = make_function(H2_symbolic, r_ij_vars)

    # Three body part need to do a bit differently
    r_i_vars = make_variables(:ri, D)
    r_j_vars = make_variables(:rj, D)
    r_k_vars = make_variables(:rk, D)

    r_ij_threebody = r_i_vars .- r_j_vars
    r_ik_threebody = r_i_vars .- r_k_vars
    r_ij_norm_threebody = sqrt(sum(x -> x^2, r_ij_threebody))
    r_ik_norm_threebody = sqrt(sum(x -> x^2, r_ik_threebody))

    #3rd order part
    pot3_symbolic = three_body_potential_nounits(pot, r_ij_threebody, r_ik_threebody, r_ij_norm_threebody, r_ik_norm_threebody)

    #Can't treat block of IFC as hessian
    H3_symbolic_ij = Matrix{FastDifferentiation.Node}(undef, D, D)
    for a in range(1,D)
        for b in range(a,D)
            H3_symbolic_ij[a,b] = derivative([pot3_symbolic], r_i_vars[a], r_j_vars[b])[1]
            H3_symbolic_ij[b,a] = H3_symbolic_ij[a,b]
        end
    end

    H3_symbolic_jk = Matrix{FastDifferentiation.Node}(undef, D, D)
    for a in range(1,D)
        for b in range(a,D)
            H3_symbolic_jk[a,b] = derivative([pot3_symbolic], r_j_vars[a], r_k_vars[b])[1]
            H3_symbolic_jk[b,a] = H3_symbolic_jk[a,b]
        end
    end

    H3_symbolic_ik = Matrix{FastDifferentiation.Node}(undef, D, D)
    for a in range(1,D)
        for b in range(a,D)
            H3_symbolic_ik[a,b] = derivative([pot3_symbolic], r_i_vars[a], r_k_vars[b])[1]
            H3_symbolic_ik[b,a] = H3_symbolic_ik[a,b]
        end
    end

    H3_exec_ij = make_function(H3_symbolic_ij, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_jk = make_function(H3_symbolic_jk, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_ik = make_function(H3_symbolic_ik, [r_i_vars; r_j_vars; r_k_vars])

    return H2_exec, H3_exec_ij, H3_exec_jk, H3_exec_ik
end