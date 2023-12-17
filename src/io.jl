export save_second_order, save_third_order,
        parse_TDEP_second_order, parse_TDEP_third_order,
        parse_ModeCode_third_order, parse_LAMMPS_dynmat, parse_LAMMPS_third_order,
        parse_xyz


# Methods to save data created with this package

function save_second_order(Φ::SecondOrderForceConstants, N_atoms, D, outpath; filename = "second_order", fmt = :JLD2)
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

function save_third_order(Ψ::ThirdOrderForceConstants, N_atoms, D, outpath; filename = "third_order", fmt = :JLD2)
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
        jldsave(filepath, F3 = Ψ.values, N_atoms = N_atoms, dimension = D, unit = string(Ψ.units))
    else
        throw(ArgumentError("Invalid format, $(fmt)"))
    end
end

function save_third_order(Ψ::SparseThirdOrder, N_atoms, D, outpath; filename = "third_order_sparse", fmt = :JLD2)

    if fmt == :JLD2
        filepath = joinpath(outpath, filename*".jld2")
        jldsave(filepath, F3 = Ψ.values, N_atoms = N_atoms, dimension = D, unit = string(Ψ.units))
    else
        throw(ArgumentError("Invalid format, $(fmt)"))
    end

end

#### Methods to load force constants from other libraries

"""
Parses the TDEP format for second order force cosntants. 

This function assumes that the force constants were calculated for the supercell or that
the unit-cell was remapped to the supercell with remap_forceconstants. This function is
not guranteed to work with the most recent version of TDEP.
"""
function parse_TDEP_second_order(ifc_path::String, N_modes, energy_units::Symbol = :REAL)

    Φ = zeros(N_modes, N_modes)

    open(ifc_path, "r") do f
    
        #First 2 lines are header
        N_atoms = parse(Int64,split(strip(readline(f)))[1])
        @assert 3*N_atoms == N_modes "Cannot handle non-monatomic case"
        r_cut = parse(Float64,split(strip(readline(f)))[1])
    
        #Next lines gives num neighbors of atom in unit cell we are looking at
        Φ_block = zeros(3,3)
        for i in 1:N_atoms
            n_neighbors, base_atom_idx_uc = parse.(Int64,split(strip(readline(f)))[[1,end-1]])
            
            for j in 1:n_neighbors
                other_atom_idx_uc = parse(Int64,split(strip(readline(f)))[1])
                lattice_vec = parse.(Float64,split(strip(readline(f)))) #Dont think this is useful?

                Φ_block[1,:] .= parse.(Float64,split(strip(readline(f)))) 
                Φ_block[2,:] .= parse.(Float64,split(strip(readline(f))))
                Φ_block[3,:] .= parse.(Float64,split(strip(readline(f))))

                Φ[3*(base_atom_idx_uc - 1) + 1 : 3*base_atom_idx_uc,
                  3*(other_atom_idx_uc - 1) + 1 : 3*other_atom_idx_uc] .= Φ_block

            end
        end
    end

    if energy_units == :REAL
        Φ .*= 23.060541945
        units = u"kcal * mol^-1 * Å^-2"
    elseif energy_units == :METAL
        units = u"eV * Å^-2"
    else
        throw(ArgumentError("Unsupported unit type: $(energy_units)"))
    end

    return DenseForceConstants(Φ, units, 0.0)

end

"""
Parses the TDEP format for third order force cosntants. 

This function assumes that the force constants were calculated for the supercell or that
the unit-cell was remapped to the supercell with remap_forceconstants. This function is
not guranteed to work with the most recent version of TDEP.
"""
function parse_TDEP_third_order(ifc_path::String, N_modes, energy_units = :REAL)
    Ψ = zeros(N_modes, N_modes, N_modes)

    open(ifc_path, "r") do f
    
        #First 2 lines are header
        N_atoms = parse(Int64,split(strip(readline(f)))[1])
        @assert 3*N_atoms == N_modes "Cannot handle non-monatomic case"
        r_cut = parse(Float64,split(strip(readline(f)))[1])
    
        #Next lines gives num neighbors of atom in unit cell we are looking at
        Ψ_block = zeros(3,3,3)
        for i in 1:N_atoms
            n_triplets, base_atom_idx_uc = parse.(Int64,split(strip(readline(f)))[[1,end-2]])
            
            for j in 1:n_triplets
                i = parse(Int64,split(strip(readline(f)))[1])
                j = parse(Int64,split(strip(readline(f)))[1])
                k = parse(Int64,split(strip(readline(f)))[1])
                lattice_vec = parse.(Float64,split(strip(readline(f)))) #Dont think this is useful?
                lattice_vec = parse.(Float64,split(strip(readline(f)))) #Dont think this is useful?
                lattice_vec = parse.(Float64,split(strip(readline(f)))) #Dont think this is useful?

                Ψ_block[1,1,:] .= parse.(Float64,split(strip(readline(f)))) 
                Ψ_block[1,2,:] .= parse.(Float64,split(strip(readline(f))))
                Ψ_block[1,3,:] .= parse.(Float64,split(strip(readline(f))))

                Ψ_block[2,1,:] .= parse.(Float64,split(strip(readline(f)))) 
                Ψ_block[2,2,:] .= parse.(Float64,split(strip(readline(f))))
                Ψ_block[2,3,:] .= parse.(Float64,split(strip(readline(f))))

                Ψ_block[3,1,:] .= parse.(Float64,split(strip(readline(f)))) 
                Ψ_block[3,2,:] .= parse.(Float64,split(strip(readline(f))))
                Ψ_block[3,3,:] .= parse.(Float64,split(strip(readline(f))))

                Ψ[3*(i - 1) + 1 : 3*i,
                  3*(j - 1) + 1 : 3*j,
                  3*(k - 1) + 1 : 3*k] .= Ψ_block

            end
        end
    end

    if energy_units == :REAL
        Ψ .*= 23.060541945
        units = u"kcal * mol^-1 * Å^-3"
    elseif energy_units == :METAL
        units = u"eV * Å^-3"
    else
        throw(ArgumentError("Unsupported unit type: $(energy_units)"))
    end

    return DenseForceConstants(Ψ, units, 0.0)
end

function parse_ModeCode_second_order(path::String, asr_path::String, N_atoms::Int, unit_system::Symbol = :REAL)
    conversion_factor = 1.0
    if unit_system == :REAL
        # 1.889 bohr = 1 Å
        # 313.75470835207074 kcal/mol = 1 Ryd
        units = u"kcal * mol^-1 * Å^-2"
        conversion_factor = (1.889725988*1.889725988)*313.75470835207074
    elseif unit_system == :METAL
        # 13.60569301 eV = 1 Ryd
        units = u"eV * Å^-2"
        conversion_factor = (1.889725988*1.889725988)*13.60569301
    else
        throw(ArgumentError("Unsupported unit_system: $(unit_system)"))
    end

    F2_file_contents = readdlm(path);
    F2 = zeros((3*N_atoms, 3*N_atoms))
    #first line is number of F3
    for line in range(2,size(F2_file_contents)[1])
        i, α, j, β, data = F2_file_contents[line,:]
        F2[Int((3*(i-1)) + α), Int((3*(j-1)) + β)] = data*conversion_factor
    end

    ASR_file_contents = readdlm(asr_path)
    for line in range(2,size(ASR_file_contents)[1])
        i, α, j, β, data = ASR_file_contents[line,:]
        F2[Int((3*(i-1)) + α), Int((3*(j-1)) + β)] = data*conversion_factor
    end

    return DenseForceConstants(F2, units, 0.0)
end

"""
Modecode assumes LAMMPS units were metal and converts to Ryd/Bohr^3.
"""
function parse_ModeCode_third_order(path::String, N_atoms::Int, unit_system::Symbol = :REAL)
    #Convert Ryd/bohr^3 to match system units
    conversion_factor = 1.0
    if unit_system == :REAL
        # 1.889 bohr = 1 Å
        # 313.75470835207074 kcal/mol = 1 Ryd
        units = u"kcal * mol^-1 * Å^-3"
        conversion_factor = (1.889725988*1.889725988*1.889725988)*313.75470835207074
    elseif unit_system == :METAL
        # 13.60569301 eV = 1 Ryd
        units = u"eV * Å^-3"
        conversion_factor = (1.889725988*1.889725988*1.889725988)*13.60569301
    else
        throw(ArgumentError("Unsupported unit_system: $(unit_system)"))
    end

    F3_file_contents = readdlm(path);

    F3 = zeros((3*N_atoms, 3*N_atoms, 3*N_atoms))
    #first line is number of F3
    for line in range(2,size(F3_file_contents)[1])
        i, α, j, β, k, γ, data = F3_file_contents[line,:]
        F3[Int((3*(i-1)) + α), Int((3*(j-1)) + β), Int((3*(k-1)) + γ)] = data*conversion_factor
    end
    return DenseForceConstants(F3, units, 0.0)
end

"""
Parses output from LAMMPS `dynamical_matrix` command in the PHONON package
"""
function parse_LAMMPS_dynmat(path::String, N_atoms::Integer, units)

    dynmat_file_contents = readdlm(path);

    dynmat = zeros((3*N_atoms, 3*N_atoms))

    row = 1

    for i in range(1,N_atoms)
        for alpha in range(1,3)
            for j in range(1,N_atoms)
                for beta in range(1,3)
                    dynmat[(3*(i-1))+alpha, (3*(j-1))+beta] = dynmat_file_contents[row,beta]
                end
                row += 1
            end
        end
    end

    if unit_system == :REAL
        units = u"kcal * g^-1 * Å^-2"
    elseif unit_system == :METAL
        units = u"eV * mol * g^-1 * Å^-2"
    else
        throw(ArgumentError("Unsupported unit_system: $(unit_system)"))
    end

    return DenseForceConstants(dynmat, units, 0.0)
end
"""
Parses output from LAMMPS `third_order` command in the PHONON package
"""
function parse_LAMMPS_third_order(path::String, N_atoms::Int, unit_system::Symbol = :REAL)

    F3_file_contents = readdlm(path);

    F3 = zeros((3*N_atoms, 3*N_atoms, 3*N_atoms))

    for line in range(1,size(F3_file_contents)[1])
        i, α, j, β, k, data... = F3_file_contents[line,:]
        for γ in range(1,3)
            F3[Int((3*(i-1)) + α), Int((3*(j-1)) + β), Int((3*(k-1)) + γ)] = data[γ]
        end
    end

    if unit_system == :REAL
        units = u"kcal * mol^-1 * Å^-3"
    elseif unit_system == :METAL
        units = u"eV * Å^-3"
    else
        throw(ArgumentError("Unsupported unit_system: $(unit_system)"))
    end

    return DenseForceConstants(F3, units, 0.0)
end


function parse_xyz(path, N_atoms, N_steps)

    file_contents = readdlm(path);
    positions = Vector{Vector{Vector{Float64}}}(undef, (N_steps,))
    ids = Vector{Vector{Int64}}(undef, (N_steps,))

    for i in 1:N_steps
        positions[i] = Vector{Float64}[]
        ids[i] = Int64[]
    end

    line_num = 0
    step_num = 1
    header_size = 2

    for line in range(1,size(file_contents)[1])
       
        if line_num == 0 || line_num == 1
            line_num += 1
            continue
        end

        id, posn... = file_contents[line,:]
        push!(positions[step_num], posn)
        push!(ids[step_num], Int64(id))

        line_num += 1

        if line_num == N_atoms + header_size
            line_num = 0
            step_num += 1
            continue
        end
    end
    return positions, ids
end