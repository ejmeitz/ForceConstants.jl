export save_second_order, save_third_order

#CHANGE THIS TO USE MULTIPLE DISPATCH 

function save_second_order(Φ::SecondOrderMatrix, N_atoms, D, outpath; filename = "second_order", fmt = :JLD2)
    @assert ((N_atoms*D) == size(Φ.values)[1] && (N_atoms*D) == size(Φ.values)[2]) "Incorrect dimensions"

    if fmt == :txt
        filepath = joinpath(outpath, filename*".txt")
        f = open(filepath, "w")
        println(f,"#$(N_atoms) atoms, Sum Φ = $(round(sum(Φ.values), sigdigits = 7)), Units: $(Φ.units[1,1])")
        for i in range(1, N_atoms)
            for j in range(1, N_atoms)
                for α in range(1,D)
                    for β in range(1,D)
                        ii = D*(i-1) + α; jj = D*(j-1) + β
                        if Φ.values[ii,jj] != 0.0
                            println(f, "$i $α $j $β $(Φ.values[ii,jj])")
                        end
                    end
                end
            end
        end

        close(f)
    elseif fmt == :JLD2
        filepath = joinpath(outpath, filename*".jld2")
        jldsave(filepath, F2 = Φ.values, N_atoms = N_atoms, dimension = D, unit = string(Φ.units))
    else
        throw(ArgumentError("Invalid format, $(fmt)"))
    end

end

function save_third_order(Ψ::ThirdOrderMatrix, N_atoms, D, outpath; filename = "third_order", fmt = :JLD2)
    @assert ((N_atoms*D) == size(Ψ.values)[1] && (N_atoms*D) == size(Ψ)[2] && (N_atoms*D) == size(Ψ.values)[3]) "Incorrect dimensions"

    if fmt == :txt
        filepath = joinpath(outpath, filename*".txt")
        f = open(filepath, "w")
        println(f,"#$(N_atoms) atoms, Sum Ψ = $(round(sum(Ψ.values),sigdigits = 7)), Units: $(Ψ.units)")
        for i in range(1, N_atoms)
            for j in range(1, N_atoms)
                for k in range(1,N_atoms)
                    for α in range(1,D)
                        for β in range(1,D)
                            for γ in range(1,D)
                                ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(k-1) + γ
                                if Ψ.values[ii,jj,kk] != 0.0
                                    println(f, "$i $α $j $β $k $γ $(Ψ.values[ii,jj,kk])")
                                end
                            end
                        end
                    end
                end
            end
        end
        close(f)
    elseif fmt == :JLD2
        filepath = joinpath(outpath, filename*".jld2")
        jldsave(filepath, F3 = Ψ, N_atoms = N_atoms, dimension = D, unit = string(Ψ.units))
    else
        throw(ArgumentError("Invalid format, $(fmt)"))
    end
end

function save_third_order(Ψ::ThirdOrderSparse, N_atoms, D, outpath; filename = "third_order_sparse", fmt = :JLD2)

    if fmt == :JLD2
        filepath = joinpath(outpath, filename*".jld2")
        jldsave(filepath, F3 = Ψ.values, N_atoms = N_atoms, dimension = D, unit = string(Ψ.units))
    else
        throw(ArgumentError("Invalid format, $(fmt)"))
    end

end