export gruneisen_parameter

"""
Calculate mode Gruneisen Parameter's from second and third order force constants
"""
function gruneisen_parameter(sys::SuperCellSystem{D}, 
        Φ::Union{SecondOrderForceConstants, Matrix{T}},
        Ψ::Union{ThirdOrderForceConstants,Array{T,3}}) where {D,T}

    N_modes = size(Φ)[1]
    N_atoms = Int64(N_modes/D)

    @assert n_atoms(sys) == N_atoms "Size of force constant matrix does not match n_atoms of system"

    dynmat = dynamical_matrix(sys, Φ)
    freqs_sq, e_vecs = get_modes(dynmat, D)

    #Pre-compute r_bar and mass normalization
    mass_norm = calculate_mass_normalization(sys)
    r_flat = ustrip.(positions_1D(sys))
    # r_bar = calculate_r_bar(sys, Φ, e_vecs, freqs_sq, mass_norm, r_flat)
    # r_bar_flat = reduce(vcat, r_bar)
    # print(r_bar_flat)
    r_bar_flat = zeros(D*N_atoms)
    r_bar_flat .+= r_flat
    
    gamma = zeros(N_modes)

    Threads.@threads for n in range(1,N_modes)
        @info "$n"
        if freqs_sq[n] != 0.0
            for i in range(1,N_atoms)
                for j in range(1,N_atoms)
                    # tmp = 0.0
                    for k in range(1,N_atoms)
                        for α in range(1,D)
                            ii = D*(i-1) + α
                            for β in range(1,D)
                                jj = D*(j-1) + β
                                for γ in range(1,D)
                                    kk = D*(k-1) + γ
                                    gamma[n] += Ψ[ii,jj,kk]*conj(e_vecs[ii,n])*e_vecs[jj,n]*r_bar_flat[kk]*mass_norm[i,j]
                                end
                            end
                        end
                    end
                    # gamma[n] += tmp*mass_norm[i,j]
                end
            end
        end
    end
    gamma ./= (-6 .*freqs_sq)

    #Remove NaNs
    gamma[freqs_sq .== 0.0] .= 0.0
    return gamma
    # nothing

end

function calculate_mass_normalization(sys::SuperCellSystem{D}) where {D}

    N_atoms = n_atoms(sys)
    mass_norm = zeros(N_atoms,N_atoms)
    Threads.@threads for i in range(1,N_atoms)
        for j in range(i,N_atoms)
            mass_norm[i,j] = 1/ustrip(sqrt(mass(sys, i)*mass(sys,j)))
            mass_norm[j,i] = mass_norm[i,j]
        end
    end

    return mass_norm
end


function calculate_r_bar(sys::SuperCellSystem{D}, Φ, e_vecs, freqs_sq, mass_norm, r_flat) where {D}

    N_atoms = n_atoms(sys)
    N_modes = size(e_vecs)[1]
    r_bar = zeros(N_atoms,D)
    
    
    Threads.@threads for i in range(1,N_atoms)
        for α in range(1,D)
            ii = D*(i-1) + α
            for n in range(1,N_modes)
                if freqs_sq[n] != 0.0
                    tmp = 0.0
                    for j in range(1,N_atoms)
                        for k in range(1,N_atoms)
                            for β in range(1,D)
                                jj = D*(j-1) + β
                                for γ in range(1,D)
                                    kk = D*(k-1) + γ
                                    tmp += Φ[jj,kk]*e_vecs[jj,n]*r_flat[kk]*mass_norm[i,j]
                                    # r_bar[i,α] += Φ[jj,kk]*e_vecs[jj,n]*r_flat[kk]*mass_norm[i,j]*e_vecs[ii,n]/freqs_sq[n]
                                end
                            end
                        end
                    end
                    if tmp > 1e-10
                        println(tmp)
                    end
                    r_bar[i,α] += tmp*conj(e_vecs[ii,n])/freqs_sq[n]
                end
            end
        end
    end

    r_bar .*= -1

    return r_bar
end


"""
Calculate mode Gruneisen Parameter's as a derivative of the frequencies.
Will accept either 2 or 4 systems and perform a central difference. For example,
with two systems, one system should have its volume increased by dV and the other decreased by dV.
The order of the input does not matter.
"""
function gruneisen_parameter(systems::Vector{<:SuperCellSystem{D}},
     pot::Potential, calc::ForceConstantCalculator, V₀, base_freqs) where {D}


    @assert length(systems) ∈ [2,4] "Must provide 2 or 4 systems for Gruneisen Parameter calculation"
    @assert allequal(n_atoms.(systems)) "Systems must have same number of atoms"

    freqs = Vector{Float64}[]

    #Sort smallest to largest volume
    # -dV, +dV for 2 systems
    # -2dV, -dV, +dV, +2dV for 4 systems
    systems = sort(systems, by = (sys -> volume(sys)))

    for sys in systems
        dynmat = dynamical_matrix(sys, pot, calc).values
        freqs_sq, _ = get_modes(dynmat, D)
        push!(freqs, sqrt.(freqs_sq))
    end

    N_modes = length(freqs[1])
    derivs = zeros(N_modes)

    if length(systems) == 2
        dV = ustrip(0.5*volume(systems[2]) - volume(systems[1]))
        Threads.@threads for n in eachindex(freqs[1])
            derivs[n] = (-freqs[1][n] + freqs[2][n])/(2*dV)
        end
    elseif length(systems) == 4
        dV = ustrip(volume(systems[4]) - volume(systems[3]))
        Threads.@threads for n in eachindex(freqs[1])
            derivs[n] = (freqs[1][n] - 8*freqs[2][n] + 8*freqs[3][n] - freqs[4][n])/(12*dV)
        end
    end

    gamma = -derivs.*(ustrip(V₀)./ustrip.(base_freqs))
    gamma[base_freqs .== 0.0] .= 0.0

    return gamma
end