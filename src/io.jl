export save_second_order, save_third_order

#CHANGE THIS TO USE MULTIPLE DISPATCH 

function save_second_order(Φ, N_atoms, D, outpath; filename = "second_order", fmt = :txt)
    @assert ((N_atoms*D) == size(Φ)[1] && (N_atoms*D) == size(Φ)[2]) "Incorrect dimensions"

    if fmt == :txt
        filepath = joinpath(outpath, filename*".txt")
        f = open(filepath, "w")
        println(f,"#$(N_atoms) atoms, Sum Φ = $(round(ustrip(sum(Φ)),sigdigits = 4)), Units: $(unit(Φ[1,1]))")
        for i in range(1, N_atoms)
            for j in range(1, N_atoms)
                for α in range(1,D)
                    for β in range(1,D)
                        ii = D*(i-1) + α; jj = D*(j-1) + β
                        if ustrip(Φ[ii,jj]) != 0.0
                            println(f, "$i $α $j $β $(ustrip(Φ[ii,jj]))")
                        end
                    end
                end
            end
        end

        close(f)
    elseif fmt == :JLD2
        throw(ArgumentError("Not implemented yet, $(fmt)"))
    else
        throw(ArgumentError("Invalid format, $(fmt)"))
    end

end

function save_third_order(Ψ, N_atoms, D, outpath; filename = "third_order", fmt = :txt)
    @assert ((N_atoms*D) == size(Ψ)[1] && (N_atoms*D) == size(Ψ)[2] && (N_atoms*D) == size(Ψ)[3]) "Incorrect dimensions"

    if fmt == :txt
        filepath = joinpath(outpath, filename*".txt")
        f = open(filepath, "w")
        println(f,"#$(N_atoms) atoms, Sum Ψ = $(round(ustrip(sum(Ψ)),sigdigits = 4)), Units: $(unit(Ψ[1,1,1]))")
        for i in range(1, N_atoms)
            for j in range(1, N_atoms)
                for k in range(1,N_atoms)
                    for α in range(1,D)
                        for β in range(1,D)
                            for γ in range(1,D)
                                ii = D*(i-1) + α; jj = D*(j-1) + β; kk = D*(k-1) + γ
                                if ustrip(Ψ[ii,jj,kk]) != 0.0
                                    println(f, "$i $α $j $β $k $γ $(ustrip(Ψ[ii,jj,kk]))")
                                end
                            end
                        end
                    end
                end
            end
        end
        close(f)
    elseif fmt == :JLD2
        throw(ArgumentError("Not implemented yet, $(fmt)"))
    else
        throw(ArgumentError("Invalid format, $(fmt)"))
    end
end


function save_third_order(Ψ::ThirdOrderSparse, N_atoms, D, outpath; filename = "third_order", fmt = :txt)


end