export gruneisen_parameter, calculate_r_bar2, calculate_r_bar3

"""
Calculate mode Gruneisen Parameter's from second and third order force constants
"""
function gruneisen_parameter(sys::SuperCellSystem{D}, 
        Φ::Union{SecondOrderForceConstants, Matrix{T}},
        Ψ::Union{ThirdOrderForceConstants,Array{T,3}}) where {D,T}

    N_modes = size(Φ)[1]
    N_atoms = Int64(N_modes/D)
    println("in gruneisen_parameter")
    @assert n_atoms(sys) == N_atoms "Size of force constant matrix does not match n_atoms of system"

    dynmat = dynamical_matrix(sys, Φ)
    freqs_sq, e_vecs = get_modes(dynmat, D)

    #Pre-compute r_bar and mass normalization
    mass_norm = calculate_mass_normalization(sys)
    r_flat = ustrip.(positions_1D(sys))
    # r_bar = calculate_r_bar(sys, Φ, e_vecs, freqs_sq, mass_norm, r_flat)

    # r_bar_flat = reduce(vcat, r_bar)
    # println("R Bar Range: $(minimum(r_bar_flat)) to $(maximum(r_bar_flat))")

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
            iα = D*(i-1) + α
            for n in range(1,N_modes)
                if freqs_sq[n] != 0.0
                    tmp = 0.0
                    for jβ in range(1,N_modes)
                        j_idx_1D = floor(Int64,(jβ-1)/D) + 1
                        for kγ in range(1,N_modes)
                            tmp += Φ[jβ,kγ]*e_vecs[jβ,n]*r_flat[kγ]*mass_norm[i,j_idx_1D]
                            # r_bar[i,α] += Φ[jβ,kγ]*e_vecs[jβ,n]*r_flat[kγ]*mass_norm[i,j_idx_1D]*e_vecs[iα,n]/freqs_sq[n]
                        end
                    end
                    r_bar[i,α] += tmp*conj(e_vecs[iα,n])/freqs_sq[n]
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

function make_grunesien_test_systems(dV, n_sys, n_uc, calc, pot)

    a0 = 5.43u"Å"
    base_crys = Diamond(a0, :Si, SVector(n_uc,n_uc,n_uc))
    base_sys = SuperCellSystem(base_crys)
    V₀ = volume(base_sys)
    da = (cbrt(V₀ + dV) - n_uc*a0)/n_uc

    dynmat = dynamical_matrix(base_sys, pot, calc).values
    freqs_sq, _ = get_modes(dynmat)
    base_freqs = sqrt.(freqs_sq)

    if n_sys == 2
        crys1 = Diamond(a0 + da, :Si, SVector(n_uc,n_uc,n_uc))
        crys2 = Diamond(a0 - da, :Si, SVector(n_uc,n_uc,n_uc))
        sys1 = SuperCellSystem(crys1)
        sys2 = SuperCellSystem(crys2)

        return [sys1, sys2], base_sys, V₀
    elseif n_sys == 4
        crys1 = Diamond(a0 + 2*da, :Si, SVector(n_uc,n_uc,n_uc))
        crys2 = Diamond(a0 + da, :Si, SVector(n_uc,n_uc,n_uc))
        crys3 = Diamond(a0 - da, :Si, SVector(n_uc,n_uc,n_uc))
        crys4 = Diamond(a0 - 2*da, :Si, SVector(n_uc,n_uc,n_uc))
        sys1 = SuperCellSystem(crys1)
        sys2 = SuperCellSystem(crys2)
        sys3 = SuperCellSystem(crys3)
        sys4 = SuperCellSystem(crys4)

        return [sys1, sys2, sys3, sys4], base_freqs, V₀
    end

end


function calculate_r_bar2(sys::SuperCellSystem{D}, Φ, e_vecs, freqs_sq,
         mass_norm, r_flat, i, α::Integer) where {D}

    N_atoms = n_atoms(sys)
    N_modes = size(e_vecs)[1]
    r_bar = zeros(N_atoms,D)
    println("HERE2")
    #Sanity checks
    @assert 3*N_atoms == N_modes

    # Threads.@threads for i in range(1,N_atoms)
    #     for α in range(1,D)
    ii = D*(i-1) + α
    for n in range(1,N_modes)
        if freqs_sq[n] != 0.0
            tmp = 0.0
            for j in range(1,N_atoms)
                for k in range(1,N_atoms)
                    for β in range(1,D)
                        jj = D*(j-1) + β;
                        for γ in range(1,D)
                            kk = D*(k-1) + γ
                            tmp += Φ[jj,kk]*e_vecs[jj,n]*r_flat[kk]*mass_norm[i,j]
                        end
                    end
                end
            end
            r_bar[i,α] += tmp*(conj(e_vecs[ii,n])/freqs_sq[n])
        end
    end
    @info "$(-r_bar[i,α])"
    #     end
    # end

    r_bar .*= -1

    return nothing
end

function calculate_r_bar3(sys::SuperCellSystem{D}, Φ, e_vecs, freqs_sq, mass_norm,
     r_flat, i, α::Integer) where {D}
    N_atoms = n_atoms(sys)
    N_modes = size(e_vecs)[1]
    r_bar = zeros(N_atoms,D)
    println("HERE")

    #Sanity checks
    @assert 3*N_atoms == N_modes
    
    ii = D*(i-1) + α
    for n in range(1,N_modes)
        if freqs_sq[n] != 0.0
            for j in range(1,N_atoms)
                for k in range(1,N_atoms)
                    for β in range(1,D)
                        for γ in range(1,D)
                            jj = D*(j-1) + β; kk = D*(k-1) + γ
                            r_bar[i,α] += Φ[jj,kk]*e_vecs[jj,n]*r_flat[kk]*mass_norm[i,j]*conj(e_vecs[ii,n])/freqs_sq[n]
                        end
                    end
                end
            end
        end
    end
    @info "$(-r_bar[i,α])"

    r_bar .*= -1

    return nothing
end