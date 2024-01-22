export dynamical_matrix, dynamical_matrix!, get_modes

### Super Cell ###

function dynamical_matrix(sys::SuperCellSystem{D}, pot::Potential, 
    calc::ForceConstantCalculator) where D

    if any(calc.r_cut .> sys.box_sizes_SC)
         @warn "Cutoff larger than L/2"
    end
    N_atoms = n_atoms(sys)

    #reuse storage from IFC2 calculation
    dynmat = second_order(sys, pot, calc)

    #Mass Weight
    Threads.@threads for i in range(1, N_atoms)
        for j in range(1, N_atoms)
            for α in range(1,D)
                for β in range(1,D)
                    ii = D*(i-1) + α
                    jj = D*(j-1) + β
                    dynmat.values[ii,jj] /=  ustrip(sqrt(mass(sys,i)*mass(sys,j)))
                end
            end
        end
    end

    #Add final units to dynamical matrix
    dynmat_unit = dynmat.units / unit(mass(sys,1))

    return DenseForceConstants(dynmat.values, dynmat_unit, 0.0)
end

function dynamical_matrix!(dynmat, sys::SuperCellSystem{D}, pot::Potential, 
    calc::ForceConstantCalculator) where D

    if any(calc.r_cut .> sys.box_sizes_SC)
        @warn "Cutoff larger than L/2"
    end
    N_atoms = n_atoms(sys)

    @assert size(dynmat) == (D*N_atoms, D*N_atoms)

    #reuse storage from IFC2 calculation
    dynmat = second_order!(dynmat, sys, pot, calc)

    #Mass Weight
    Threads.@threads for i in range(1, N_atoms)
        for j in range(1, N_atoms)
            for α in range(1,D)
                for β in range(1,D)
                    ii = D*(i-1) + α
                    jj = D*(j-1) + β
                    dynmat[ii,jj] /=  ustrip(sqrt(mass(sys,i)*mass(sys,j)))
                end
            end
        end
    end

    return dynmat
end

"""
Builds new dynamical matrix from existing second order force constants
"""
function dynamical_matrix(sys::SuperCellSystem{D}, Φ::SecondOrderForceConstants) where {D}

    dynmat = similar(Φ)
    dynmat .= Φ.values

    N_atoms = n_atoms(sys)

    #Mass Weight
    Threads.@threads for i in range(1, N_atoms)
        for j in range(1, N_atoms)
            for α in range(1,D)
                for β in range(1,D)
                    ii = D*(i-1) + α
                    jj = D*(j-1) + β
                    dynmat[ii,jj] /=  ustrip(sqrt(mass(sys,i)*mass(sys,j)))
                end
            end
        end
    end

    #Add final units to dynamical matrix
    dynmat_unit = Φ.units / unit(mass(sys,1))

    return DenseForceConstants(dynmat, dynmat_unit, 0.0)
end

function dynamical_matrix(sys::SuperCellSystem{D}, Φ::Matrix) where {D}
    N_atoms = n_atoms(sys)
    dynmat = similar(Φ)

    #Mass Weight
    Threads.@threads for i in range(1, N_atoms)
        for j in range(i, N_atoms)
            for α in range(1,D)
                for β in range(1,D)
                    ii = D*(i-1) + α
                    jj = D*(j-1) + β
                    dynmat[ii,jj] =  Φ[ii,jj]/ustrip(sqrt(mass(sys,i)*mass(sys,j)))
                    dynmat[jj,ii] =  dynmat[ii,jj]
                end
            end
        end
    end
    
    return dynmat

end

### Get Modes ###
"""
get_modes(dynmat::SecondOrderMatrix, num_rigid_translation = 3)

This function automatically zeros out the 3 modes whose frequency is 
nearest to zero which are assumed to be the rigid translation modes.
To change this pass the kwarg `num_rigid_translation`.
"""
function get_modes(dynmat::SecondOrderForceConstants, num_rigid_translation = 3)

    eig_stuff = eigen(Hermitian(dynmat.values))
    freqs_sq = eig_stuff.values
    idx_rt = sortperm(abs.(freqs_sq))
    freqs_sq[idx_rt[1:num_rigid_translation]] .= 0.0

    return freqs_sq, eig_stuff.vectors
end

function get_modes(dynmat::Matrix, num_rigid_translation = 3)

    eig_stuff = eigen(Hermitian(dynmat))
    freqs_sq = eig_stuff.values
    idx_rt = sortperm(abs.(freqs_sq))
    freqs_sq[idx_rt[1:num_rigid_translation]] .= 0.0

    return freqs_sq, eig_stuff.vectors
end

