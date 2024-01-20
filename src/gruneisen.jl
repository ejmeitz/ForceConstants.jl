export gruneisen_parameter

function gruneisen_parameter(sys::SuperCellSystem{D}, 
        Φ::SecondOrderForceConstants, Ψ::ThirdOrderForceConstants) where {D}

    N_modes = size(Φ.values)[1]
    N_atoms = Int64(N_modes/D)

    @assert n_atoms(sys) == N_atoms "Size of force constant matrix does not match n_atoms of system"

    dynmat = dynamical_matrix(sys, Φ)
    freqs_sq, e_vecs = get_modes(dynmat, D)

    #Pre-compute r_bar and mass normalization
    mass_norm = calculate_mass_normalization(sys)
    r_flat = ustrip.(positions_1D(sys))
    r_bar = calculate_r_bar(sys, Φ, e_vecs, freqs_sq, mass_norm, r_flat)
    r_bar_flat = reduce(vcat, r_bar)

    r_bar_flat .+= r_flat
    
    gamma = zeros(N_modes)

    Threads.@threads for n in range(1,N_modes)
        @info "$n"
        for i in range(1,N_atoms)
            for j in range(1,N_atoms)
                for k in range(1,N_atoms)
                    for α in range(1,D)
                        ii = D*(i-1) + α
                        for β in range(1,D)
                            jj = D*(j-1) + β
                            for γ in range(1,D)
                                kk = D*(k-1) + γ
                                gamma[n] += Ψ[ii,jj,kk]*e_vecs[ii,n]*e_vecs[jj,n]*r_bar_flat[kk]
                            end
                        end
                    end
                end
                gamma[n] *= mass_norm[i,j]
            end
        end
    end
    gamma *= -1/(6 .*freqs_sq)

    return gamma

end

function calculate_mass_normalization(sys::SuperCellSystem{D}) where {D}

    N_atoms = n_atoms(sys)
    mass_norm = zeros(N_atoms,N_atoms)
    Threads.@threads for i in range(1,N_atoms)
        for j in range(i,N_atoms)
            mass_norm[i,j] = 1/ustrip(sqrt(mass(sys, i)*mass(sys,j)))
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
            for n in range(1,N_modes)
                if freqs_sq[n] != 0.0
                    for j in range(1,N_atoms)
                        for k in range(1,N_atoms)
                            for β in range(1,D)
                                jj = D*(j-1) + β
                                for γ in range(1,D)
                                    kk = D*(k-1) + γ
                                    r_bar[i,α] += Φ[jj,kk]*e_vecs[jj,n]*r_flat[kk]
                                end
                            end
                        end
                        r_bar[i,α] *= mass_norm[i,j]
                    end
                    r_bar[i,α] *= e_vecs[D*(i-1) + α,n]/freqs_sq[n]
                end
            end
        end
    end

    r_bar .*= -1

    return r_bar
end



