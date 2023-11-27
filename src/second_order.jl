export dynamicalMatrix, second_order_IFC, get_modes

### Super Cell ###

function dynamicalMatrix(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    @assert all(pot.r_cut .< sys.box_sizes_SC) "Cutoff larger than L/2"
    N_atoms = n_atoms(sys)

    #reuse storage from IFC2 calculation
    dynmat = second_order_IFC(sys, pot, tol)

    #Mass Weight
    for i in range(1, N_atoms)
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

    return DenseForceConstants(dynmat.values, dynmat_unit, tol)
end

function second_order_IFC(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)

    #Loop block matricies above diagonal (not including diagonal)
    for i in range(1,N_atoms)
        for j in range(i+1,N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)
            dist = norm(rᵢⱼ)

            if dist < pot.r_cut
                for α in range(1,D)
                    for β in range(1,D)

                        #Calculate IFC2 index
                        ii = D*(i-1) + α
                        jj = D*(j-1) + β

                        IFC2[ii,jj] = -ustrip(ϕ₂(pot, dist, rᵢⱼ, α, β))
                        IFC2[jj,ii] = IFC2[ii,jj]
                    end
                end
            end
        end
    end

    #Acoustic Sum Rule -- Loop D x D block matricies on diagonal of IFC2
    for i in range(1, N_atoms) # index of block matrix
        for α in range(1,D)
            for β in range(1,D)
                ii = D*(i-1) + α
                jj = D*(i-1) + β # i == j because we're on diagonal
                IFC2[ii,jj] = -1*sum(IFC2[ii, β:D:end])
            end
        end
    end

    IFC2 = apply_tols!(IFC2,tol)

    return DenseForceConstants(IFC2, energy_unit(pot) / length_unit(pot)^2, tol)

end

function ϕ₂(pot::PairPotential, r_norm, r_jk_j′k′, α, β)
    rᵢⱼ = r_jk_j′k′

    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    return (α == β) ? ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm))) + (Φ′/r_norm) :
                    ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm)))
end

"""
Calculates the second order force constant matrix for a three-body potential. This method
uses automatic differentiation to calculate the force constants and should be checked with
the provided `check_ifc` or manually using the finite difference methods in `finite_diff_ifc.jl`.
"""
function second_order_IFC(sys::SuperCellSystem{D}, pot::ThreeBodyPotential, tol) where D
    vars = make_variables(:r, D)
    r_norm = sqrt(sum(x -> x^2, vars))
    pot_symbolic = potential_nounits(pot, r_norm)
    H_symbolic = hessian(pot_symbolic, vars)
    H_exec = make_function(H_symbolic, vars)

    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)

    Threads.@threads for i in range(1, N_atoms)
        for j in range(i + 1, N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror!(rᵢⱼ, sys.box_sizes_SC)

            ij_block = H_exec(ustrip.(rᵢⱼ))

            IFC2[D*(i-1) + 1 : D*(i-1) + D, D*(j-1) + 1 : D*(j-1) + D] .= ij_block
            IFC2[D*(j-1) + 1 : D*(j-1) + D, D*(i-1) + 1 : D*(i-1) + D] .= ij_block

        end
    end

    #Acoustic Sum Rule
    Threads.@threads for i in range(1, N_atoms) # index of block matrix
        for α in range(1,D)
            for β in range(1,D)
                ii = D*(i-1) + α
                jj = D*(i-1) + β # i == j because we're on diagonal
                IFC2[ii,jj] = -1*sum(IFC2[ii, β:D:end])
            end
        end
    end

    IFC2 = apply_tols!(IFC2,tol)

    return SecondOrderMatrix(IFC2, energy_unit(pot) / length_unit(pot)^2, tol)
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

