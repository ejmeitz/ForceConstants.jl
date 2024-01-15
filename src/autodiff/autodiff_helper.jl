function three_body_second_derivs(pot::StillingerWeberSilicon, D)

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

    H3_symbolic_ik = Matrix{FastDifferentiation.Node}(undef, D, D)
    for a in range(1,D)
        for b in range(a,D)
            H3_symbolic_ik[a,b] = derivative([pot3_symbolic], r_i_vars[a], r_k_vars[b])[1]
            H3_symbolic_ik[b,a] = H3_symbolic_ik[a,b]
        end
    end

    # H3_symbolic_jk = Matrix{FastDifferentiation.Node}(undef, D, D)
    # for a in range(1,D)
    #     for b in range(a,D)
    #         H3_symbolic_jk[a,b] = derivative([pot3_symbolic], r_j_vars[a], r_k_vars[b])[1]
    #         H3_symbolic_jk[b,a] = H3_symbolic_jk[a,b]
    #     end
    # end

    H3_exec_ij = make_function(H3_symbolic_ij, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_ik = make_function(H3_symbolic_ik, [r_i_vars; r_j_vars; r_k_vars])
    #H3_exec_jk = make_function(H3_symbolic_jk, [r_i_vars; r_j_vars; r_k_vars])

    return H2_exec, H3_exec_ij, H3_exec_ik#, H3_exec_jk
end

function three_body_third_derivs(pot::StillingerWeberSilicon, D)

    #2nd order part
    r_ij_vars = make_variables(:rij, D) #rᵢ - rⱼ
    r_ij_norm = sqrt(sum(x -> x^2, r_ij_vars))

    pot2_symbolic = pair_potential_nounits(pot, r_ij_norm)
    H_symbolic = hessian(pot2_symbolic, r_ij_vars)

    third_order_symbolic = reshape(jacobian(vec(H_symbolic), r_ij_vars),(D,D,D))
    H2_exec = make_function(third_order_symbolic, r_ij_vars)

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
    H3_symbolic_iij = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_iik = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_ijj = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_ijk = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_ikk = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_jjk = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_jkk = Array{FastDifferentiation.Node}(undef, D, D, D)
       

    for a in range(1,D)
        for b in range(a,D)
            for c in range(b,D)
                H3_symbolic_iij[a,b,c] = derivative([pot3_symbolic], r_i_vars[a], r_i_vars[b], r_j_vars[c])[1]
                H3_symbolic_iij[a,c,b] = H3_symbolic_iij[a,b,c]
                H3_symbolic_iij[b,c,a] = H3_symbolic_iij[a,b,c]
                H3_symbolic_iij[b,a,c] = H3_symbolic_iij[a,b,c]
                H3_symbolic_iij[c,a,b] = H3_symbolic_iij[a,b,c]
                H3_symbolic_iij[c,b,a] = H3_symbolic_iij[a,b,c]
            end
        end
    end

    for a in range(1,D)
        for b in range(a,D)
            for c in range(b,D)
                H3_symbolic_iik[a,b,c] = derivative([pot3_symbolic], r_i_vars[a], r_i_vars[b], r_k_vars[c])[1]
                H3_symbolic_iik[a,c,b] = H3_symbolic_iik[a,b,c]
                H3_symbolic_iik[b,c,a] = H3_symbolic_iik[a,b,c]
                H3_symbolic_iik[b,a,c] = H3_symbolic_iik[a,b,c]
                H3_symbolic_iik[c,a,b] = H3_symbolic_iik[a,b,c]
                H3_symbolic_iik[c,b,a] = H3_symbolic_iik[a,b,c]
            end
        end
    end

    for a in range(1,D)
        for b in range(a,D)
            for c in range(b,D)
                H3_symbolic_ijj[a,b,c] = derivative([pot3_symbolic], r_i_vars[a], r_j_vars[b], r_j_vars[c])[1]
                H3_symbolic_ijj[a,c,b] = H3_symbolic_ijj[a,b,c]
                H3_symbolic_ijj[b,c,a] = H3_symbolic_ijj[a,b,c]
                H3_symbolic_ijj[b,a,c] = H3_symbolic_ijj[a,b,c]
                H3_symbolic_ijj[c,a,b] = H3_symbolic_ijj[a,b,c]
                H3_symbolic_ijj[c,b,a] = H3_symbolic_ijj[a,b,c]
            end
        end
    end

    for a in range(1,D)
        for b in range(a,D)
            for c in range(b,D)
                H3_symbolic_ijk[a,b,c] = derivative([pot3_symbolic], r_i_vars[a], r_j_vars[b], r_k_vars[c])[1]
                H3_symbolic_ijk[a,c,b] = H3_symbolic_ijk[a,b,c]
                H3_symbolic_ijk[b,c,a] = H3_symbolic_ijk[a,b,c]
                H3_symbolic_ijk[b,a,c] = H3_symbolic_ijk[a,b,c]
                H3_symbolic_ijk[c,a,b] = H3_symbolic_ijk[a,b,c]
                H3_symbolic_ijk[c,b,a] = H3_symbolic_ijk[a,b,c]
            end
        end
    end

    for a in range(1,D)
        for b in range(a,D)
            for c in range(b,D)
                H3_symbolic_ikk[a,b,c] = derivative([pot3_symbolic], r_i_vars[a], r_k_vars[b], r_k_vars[c])[1]
                H3_symbolic_ikk[a,c,b] = H3_symbolic_ikk[a,b,c]
                H3_symbolic_ikk[b,c,a] = H3_symbolic_ikk[a,b,c]
                H3_symbolic_ikk[b,a,c] = H3_symbolic_ikk[a,b,c]
                H3_symbolic_ikk[c,a,b] = H3_symbolic_ikk[a,b,c]
                H3_symbolic_ikk[c,b,a] = H3_symbolic_ikk[a,b,c]
            end
        end
    end

    for a in range(1,D)
        for b in range(a,D)
            for c in range(b,D)
                H3_symbolic_jjk[a,b,c] = derivative([pot3_symbolic], r_j_vars[a], r_j_vars[b], r_k_vars[c])[1]
                H3_symbolic_jjk[a,c,b] = H3_symbolic_jjk[a,b,c]
                H3_symbolic_jjk[b,c,a] = H3_symbolic_jjk[a,b,c]
                H3_symbolic_jjk[b,a,c] = H3_symbolic_jjk[a,b,c]
                H3_symbolic_jjk[c,a,b] = H3_symbolic_jjk[a,b,c]
                H3_symbolic_jjk[c,b,a] = H3_symbolic_jjk[a,b,c]
            end
        end
    end

    for a in range(1,D)
        for b in range(a,D)
            for c in range(b,D)
                H3_symbolic_jkk[a,b,c] = derivative([pot3_symbolic], r_j_vars[a], r_k_vars[b], r_k_vars[c])[1]
                H3_symbolic_jkk[a,c,b] = H3_symbolic_jkk[a,b,c]
                H3_symbolic_jkk[b,c,a] = H3_symbolic_jkk[a,b,c]
                H3_symbolic_jkk[b,a,c] = H3_symbolic_jkk[a,b,c]
                H3_symbolic_jkk[c,a,b] = H3_symbolic_jkk[a,b,c]
                H3_symbolic_jkk[c,b,a] = H3_symbolic_jkk[a,b,c]
            end
        end
    end


    H3_exec_iij = make_function(H3_symbolic_iij, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_iik = make_function(H3_symbolic_iik, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_ijj = make_function(H3_symbolic_ijj, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_ijk = make_function(H3_symbolic_ijk, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_ikk = make_function(H3_symbolic_ikk, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_jjk = make_function(H3_symbolic_jjk, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_jkk = make_function(H3_symbolic_jkk, [r_i_vars; r_j_vars; r_k_vars])

    return H2_exec, H3_exec_iij, H3_exec_iik,
             H3_exec_ijj, H3_exec_ijk, H3_exec_ikk, H3_exec_jjk, H3_exec_jkk
end

function three_body_third_derivs_symbolic(pot::StillingerWeberSilicon, D)

    @variables ri[1:D] rj[1:D] rk[1:D]

    r_ij = ri - rj
    r_ik = ri - rk

    dist_ij = norm(r_ij)
    dist_ik = norm(r_ik)


    pot2_symbolic = pair_potential_nounits(pot, dist_ij)
    H_symbolic = hessian(pot2_symbolic, r_ij)

    H2_exec = reshape(jacobian(H_symbolic, r_ij),(D,D,D))


end