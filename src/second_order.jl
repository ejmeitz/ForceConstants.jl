export dynamicalMatrix, second_order_IFC, get_modes, second_order_finite_diff

struct SecondOrderMatrix{V,U,T}
    values::Array{V,2}
    units::U
    tol::T
end


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

    return SecondOrderMatrix(dynmat.values, dynmat_unit, tol)
end

function second_order_IFC(sys::SuperCellSystem{D}, pot::PairPotential, tol) where D
    N_atoms = n_atoms(sys)
    IFC2 = zeros(D*N_atoms,D*N_atoms)

    #Loop block matricies above diagonal (not including diagonal)
    for i in range(1,N_atoms)
        for j in range(i+1,N_atoms)

            rᵢⱼ = sys.atoms.position[i] .- sys.atoms.position[j]
            rᵢⱼ = nearest_mirror(rᵢⱼ, sys.box_sizes_SC)
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

    return SecondOrderMatrix(IFC2, unit(pot.ϵ / pot.σ^2), tol)

end

### Utility Functions ###

function ϕ₂(pot::PairPotential, r_norm, r_jk_j′k′, α, β)
    rᵢⱼ = r_jk_j′k′

    Φ′ = potential_first_deriv(pot, r_norm)
    Φ′′ = potential_second_deriv(pot, r_norm)
    return (α == β) ? ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm))) + (Φ′/r_norm) :
                    ((rᵢⱼ[α]*rᵢⱼ[β]/(r_norm^2))*(Φ′′ - (Φ′/r_norm)))
end

### Get Modes ###
"""
get_modes(dynmat::SecondOrderMatrix, num_rigid_translation = 3)

This function automatically zeros out the 3 modes whose frequency is 
nearest to zero which are assumed to be the rigid translation modes.
To change this pass the kwarg `num_rigid_translation`.
"""
function get_modes(dynmat::SecondOrderMatrix, num_rigid_translation = 3)

    eig_stuff = eigen(Hermitian(dynmat.values))
    freqs_sq = eig_stuff.values
    idx_rt = sortperm(abs.(freqs_sq))
    freqs_sq[idx_rt[1:num_rigid_translation]] .= 0.0

    return freqs_sq, eig_stuff.vectors
end


function second_order_finite_diff(sys_eq::SuperCellSystem{3}, pot::PairPotential, atom_idxs, cartesian_idxs;
    r_cut = pot.r_cut, h = 0.04*0.5291772109u"Å")

   h = uconvert(unit(pot.σ), h)
   N_atoms = n_atoms(sys_eq)


   @assert length(atom_idxs) == 2
   @assert length(cartesian_idxs) == 2
   @assert all(atom_idxs .<= N_atoms) && all(atom_idxs .>= 1) "Atom indexes out of range, must be in 1:$(N_atoms)"
   @assert all(cartesian_idxs .<= 3) && all(cartesian_idxs .>= 1) "Cartesian indices must be 1, 2, or 3"

   energy_unit = zero(potential(pot,1u"Å")) 

   #Make mutable #& change SimpleCrystals to not use SVector
   posns = [Vector(a) for a in positions(sys_eq)]

   if atom_idxs[1] != atom_idxs[2]

        energies = zeros(4)*energy_unit
        combos = [[h,h],[-h,-h],[h,-h],[-h,h]]

        for (c,combo) in enumerate(combos)

            posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]
            posns[atom_idxs[2]][cartesian_idxs[2]] += combo[2]

            energies[c] = energy_loop(pot, posns, energy_unit, r_cut, sys_eq.box_sizes_SC, N_atoms)

            #Un-modify
            posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
            posns[atom_idxs[2]][cartesian_idxs[2]] -= combo[2]
        end

        return (1/(4*h^2))*(energies[1] + energies[2] - energies[3] - energies[4])
    else
        energies = zeros(3)*energy_unit
        combos = [h,zero(pot.σ),-h]

        for (c,combo) in enumerate(combos)
            posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]

            energies[c] = energy_loop(pot, posns, energy_unit, r_cut, sys_eq.box_sizes_SC, N_atoms)

            posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
        end

        return (1/(h^2))*(energies[1] - 2*energies[2] + energies[3])

    end
end