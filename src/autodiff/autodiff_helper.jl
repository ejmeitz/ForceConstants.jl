export three_body_third_nodes, three_body_third_derivs, three_body_second_derivs

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
        for b in range(1,D)
            H3_symbolic_ij[a,b] = derivative([pot3_symbolic], r_i_vars[a], r_j_vars[b])[1]
            # H3_symbolic_ij[b,a] = H3_symbolic_ij[a,b]
        end
    end

    H3_symbolic_ik = Matrix{FastDifferentiation.Node}(undef, D, D)
    for a in range(1,D)
        for b in range(1,D)
            H3_symbolic_ik[a,b] = derivative([pot3_symbolic], r_i_vars[a], r_k_vars[b])[1]
            # H3_symbolic_ik[b,a] = H3_symbolic_ik[a,b]
        end
    end

    #*Cant assume symmetry on this one, but can for others?
    H3_symbolic_jk = Matrix{FastDifferentiation.Node}(undef, D, D)
    for a in range(1,D)
        for b in range(1,D)
            H3_symbolic_jk[a,b] = derivative([pot3_symbolic], r_j_vars[a], r_k_vars[b])[1]
            # H3_symbolic_jk[b,a] = H3_symbolic_jk[a,b]
        end
    end

    H3_exec_ij = make_function(H3_symbolic_ij, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_ik = make_function(H3_symbolic_ik, [r_i_vars; r_j_vars; r_k_vars])
    H3_exec_jk = make_function(H3_symbolic_jk, [r_i_vars; r_j_vars; r_k_vars])

    return H2_exec, H3_exec_ij, H3_exec_ik, H3_exec_jk
end

function three_body_third_nodes(pot::StillingerWeberSilicon, D)

    #2nd order part
    r_ij_vars = make_variables(:rij, D) #rᵢ - rⱼ
    r_ij_norm = sqrt(sum(x -> x^2, r_ij_vars))

    pot2_symbolic = pair_potential_nounits(pot, r_ij_norm)
    H2_symbolic_ij = hessian(pot2_symbolic, r_ij_vars)

    H3_symbolic_ij = reshape(jacobian(vec(H2_symbolic_ij), r_ij_vars),(D,D,D))
    # H3_exec_ij = make_function(H3_symbolic_ij, r_ij_vars)

    # Three body part need to do a bit differently
    r_i_vars = make_variables(:ri, D)
    r_j_vars = make_variables(:rj, D)
    r_k_vars = make_variables(:rk, D)

    r_ij_threebody = r_i_vars .- r_j_vars
    r_ik_threebody = r_i_vars .- r_k_vars
    r_ij_norm_threebody = sqrt(sum(x -> x^2, r_ij_threebody))
    r_ik_norm_threebody = sqrt(sum(x -> x^2, r_ik_threebody))

    #Can't treat block of IFC as hessian
    H3_symbolic_iij = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_iik = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_ijj = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_ijk = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_ikk = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_jjk = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_jkk = Array{FastDifferentiation.Node}(undef, D, D, D)
       
    #3rd order part
    pot3_symbolic = three_body_potential_nounits(pot, r_ij_threebody, r_ik_threebody, r_ij_norm_threebody, r_ik_norm_threebody)

    for a in range(1,D)
        for b in range(1,D)
            for c in range(1,D)
                # H3_symbolic_iij[a,b,c] = derivative([pot3_symbolic], r_i_vars[a], r_i_vars[b], r_j_vars[c])[1]
                # H3_symbolic_iik[a,b,c] = derivative([pot3_symbolic], r_i_vars[a], r_i_vars[b], r_k_vars[c])[1]
                # H3_symbolic_ijj[a,b,c] = derivative([pot3_symbolic], r_i_vars[a], r_j_vars[b], r_j_vars[c])[1]
                # H3_symbolic_ijk[a,b,c] = derivative([pot3_symbolic], r_i_vars[a], r_j_vars[b], r_k_vars[c])[1]
                # H3_symbolic_ikk[a,b,c] = derivative([pot3_symbolic], r_i_vars[a], r_k_vars[b], r_k_vars[c])[1]
                # H3_symbolic_jjk[a,b,c] = derivative([pot3_symbolic], r_j_vars[a], r_j_vars[b], r_k_vars[c])[1]
                # H3_symbolic_jkk[a,b,c] = derivative([pot3_symbolic], r_j_vars[a], r_k_vars[b], r_k_vars[c])[1]

                
                H3_symbolic_iij[a,b,c] = derivative([pot3_symbolic], r_j_vars[c], r_i_vars[b], r_i_vars[a])[1]
                H3_symbolic_iik[a,b,c] = derivative([pot3_symbolic], r_k_vars[c], r_i_vars[b], r_i_vars[a])[1]
                H3_symbolic_ijj[a,b,c] = derivative([pot3_symbolic], r_j_vars[c], r_j_vars[b], r_i_vars[a])[1]
                H3_symbolic_ijk[a,b,c] = derivative([pot3_symbolic], r_k_vars[c], r_j_vars[b], r_i_vars[a])[1]
                H3_symbolic_ikk[a,b,c] = derivative([pot3_symbolic], r_k_vars[c], r_k_vars[b], r_i_vars[a])[1]
                H3_symbolic_jjk[a,b,c] = derivative([pot3_symbolic], r_k_vars[c], r_j_vars[b], r_j_vars[a])[1]
                H3_symbolic_jkk[a,b,c] = derivative([pot3_symbolic], r_k_vars[c], r_k_vars[b], r_j_vars[a])[1]

                # H3_symbolic_iij[a,b,c] = derivative([pot3_symbolic], r_j_vars[c], r_i_vars[a], r_i_vars[b])[1]
                # H3_symbolic_iik[a,b,c] = derivative([pot3_symbolic], r_k_vars[c], r_i_vars[a], r_i_vars[b])[1]
                # H3_symbolic_ijj[a,b,c] = derivative([pot3_symbolic], r_j_vars[c], r_j_vars[b], r_i_vars[a])[1]
                # H3_symbolic_ijk[a,b,c] = derivative([pot3_symbolic], r_k_vars[c], r_j_vars[b], r_i_vars[a])[1]
                # H3_symbolic_ikk[a,b,c] = derivative([pot3_symbolic], r_k_vars[c], r_k_vars[b], r_i_vars[a])[1]
                # H3_symbolic_jjk[a,b,c] = derivative([pot3_symbolic], r_k_vars[c], r_j_vars[a], r_j_vars[b])[1]
                # H3_symbolic_jkk[a,b,c] = derivative([pot3_symbolic], r_k_vars[c], r_k_vars[b], r_j_vars[a])[1]
            end
        end
    end

    r_vars = [r_i_vars; r_j_vars; r_k_vars]

    return Dict("ij" => (H3_symbolic_ij, r_ij_vars) , "iij" => (H3_symbolic_iij, r_vars), "iik" => (H3_symbolic_iik, r_vars),
                "ijj" => (H3_symbolic_ijj, r_vars), "ijk" => (H3_symbolic_ijk, r_vars), "ikk" => (H3_symbolic_ikk, r_vars),
                "jjk" => (H3_symbolic_jjk, r_vars), "jkk" => (H3_symbolic_jkk, r_vars))

    # H3_exec_iij = make_function(H3_symbolic_iij, [r_i_vars; r_j_vars; r_k_vars])
    # H3_exec_iik = make_function(H3_symbolic_iik, [r_i_vars; r_j_vars; r_k_vars])
    # H3_exec_ijj = make_function(H3_symbolic_ijj, [r_i_vars; r_j_vars; r_k_vars])
    # H3_exec_ijk = make_function(H3_symbolic_ijk, [r_i_vars; r_j_vars; r_k_vars])
    # H3_exec_ikk = make_function(H3_symbolic_ikk, [r_i_vars; r_j_vars; r_k_vars])
    # H3_exec_jjk = make_function(H3_symbolic_jjk, [r_i_vars; r_j_vars; r_k_vars])
    # H3_exec_jkk = make_function(H3_symbolic_jkk, [r_i_vars; r_j_vars; r_k_vars])

    # return H2_exec, H3_exec_iij, H3_exec_iik,
    #          H3_exec_ijj, H3_exec_ijk, H3_exec_ikk, H3_exec_jjk, H3_exec_jkk
end

function three_body_third_derivs(pot::StillingerWeberSilicon, D)
    deriv_nodes = three_body_third_nodes(pot, D)
    return Dict(key => make_function(dn...) for (key,dn) in deriv_nodes)
end

function save_derivative(f::Array{FastDifferentiation.Node},
     vars::Vector{FastDifferentiation.Node}, base_path::String,
     filename::String)

    expr = FastDifferentiation.make_Expr(f, vars, false, true)
    @assert endswith(filename, ".jl") "Incorrect file extension, must be .jl"
    
    write(joinpath(base_path, filename), string(expr))

end


function save_SW_third_derivs(pot::StillingerWeberSilicon, D, base_path::String)

    deriv_nodes = three_body_third_nodes(pot, D)

    for (key, dn) in deriv_nodes
        filename = "H3_exec_$(key).jl"
        save_derivative(dn[1], dn[2], base_path, filename)
    end
end


function deriv_check(pot::StillingerWeberSilicon, D)
    
    # Three body part need to do a bit differently
    r_i_vars = make_variables(:ri, D)
    r_j_vars = make_variables(:rj, D)
    r_k_vars = make_variables(:rk, D)

    r_ij_threebody = r_i_vars .- r_j_vars
    r_ik_threebody = r_i_vars .- r_k_vars
    r_ij_norm_threebody = sqrt(sum(x -> x^2, r_ij_threebody))
    r_ik_norm_threebody = sqrt(sum(x -> x^2, r_ik_threebody))

    #Can't treat block of IFC as hessian
    H3_symbolic_iij = Array{FastDifferentiation.Node}(undef, D, D, D)
    H3_symbolic_iji = Array{FastDifferentiation.Node}(undef, D, D, D)
       
    #3rd order part
    pot3_symbolic = three_body_potential_nounits(pot, r_ij_threebody, r_ik_threebody, r_ij_norm_threebody, r_ik_norm_threebody)

    for a in range(1,D)
        for b in range(1,D)
            for c in range(1,D)
                H3_symbolic_iij[a,b,c] = derivative([pot3_symbolic], r_j_vars[c], r_i_vars[b], r_i_vars[a])[1]
                H3_symbolic_iji[a,b,c] = derivative([pot3_symbolic], r_i_vars[c], r_j_vars[b], r_i_vars[a])[1]
            end
        end
    end
    r_vars = [r_i_vars; r_j_vars; r_k_vars]
    H3_exec_iij = make_function(H3_symbolic_iij, r_vars)
    H3_exec_iji = make_function(H3_symbolic_iji, r_vars)

    return H3_exec_iij, H3_exec_iji

end


# function MWE(r_test)

#     r_i_vars = make_variables(:ri, D)
#     r_j_vars = make_variables(:rj, D)
#     r_k_vars = make_variables(:rk, D)

#     H3_symbolic_iij = Array{FastDifferentiation.Node}(undef, D, D, D)
#     H3_symbolic_iji = Array{FastDifferentiation.Node}(undef, D, D, D)

#     function my_norm(r_arr)
#         return sqrt(sum(x -> x^2, r_arr))
#     end

#     function f(r_arr)
#         return my_norm(r_arr[1:3]) + my_norm(r_arr[4:6]) + my_norm(r_arr[7:9])
#     end

#     f_symbolic = f([r_i_vars; r_j_vars; r_k_vars])

#     for a in range(1,D)
#         for b in range(1,D)
#             for c in range(1,D)
#                 H3_symbolic_iij[a,b,c] = derivative([f_symbolic], r_j_vars[c], r_i_vars[b], r_i_vars[a])[1]
#                 H3_symbolic_iji[a,b,c] = derivative([f_symbolic], r_i_vars[c], r_j_vars[b], r_i_vars[a])[1]
#             end
#         end
#     end

#     r_vars = [r_i_vars; r_j_vars; r_k_vars]
#     H3_exec_iij = make_function(H3_symbolic_iij, r_vars)
#     H3_exec_iji = make_function(H3_symbolic_iji, r_vars)

#     block_iij = H3_exec_iij(r_test)
#     block_iji = H3_exec_iji(r_test)

#     return block_iij == permutedims(block_iji, (1,3,2))

# end