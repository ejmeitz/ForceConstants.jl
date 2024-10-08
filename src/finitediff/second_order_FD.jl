export second_order_finite_diff_single

function second_order!(IFC2::Matrix{T}, sys_eq::SuperCellSystem{D}, pot::Potential,
    calc::FiniteDiffCalculator; n_threads = Threads.nthreads()) where {T,D}

    N_atoms = n_atoms(sys_eq)

    r_cut_sq = calc.r_cut*calc.r_cut

    @tasks for i in range(1, N_atoms)
        @set ntasks = n_threads
        @local rᵢⱼ = similar(sys_eq.atoms.position[1])
        #Give each thread a unique copy of positions
        #* has to be a more efficient way to make this thread safe
        posns = deepcopy(positions(sys_eq))
        for j in range(i+1, N_atoms)

            rᵢⱼ .= sys_eq.atoms.position[i] .- sys_eq.atoms.position[j]
            nearest_mirror!(rᵢⱼ, sys_eq.box_sizes_SC)
            dist_ij_sq = sum(x -> x^2, rᵢⱼ)
            # dist_ij_sq = 0.0u"Å^2"

            if dist_ij_sq < r_cut_sq
                for α in range(1,D)
                    ii = D*(i-1) + α
                    for β in range(1,D)
                        jj = D*(j-1) + β
                        IFC2[ii,jj] = ustrip(second_order_finite_diff_single(sys_eq, pot, [i,j], [α,β], calc.h, posns)) #* remove alloc [i,j]
                        IFC2[jj,ii] = IFC2[ii,jj]
                    end
                end
            end

        end 
    end

    ASR!(IFC2, N_atoms, D; n_threads = n_threads)

    return IFC2

end

function second_order_finite_diff_single(sys_eq::SuperCellSystem{3}, pot::Potential, atom_idxs,
     cartesian_idxs, h, posns)

   h = uconvert(length_unit(pot), h)
   N_atoms = n_atoms(sys_eq)


#    @assert length(atom_idxs) == 2
#    @assert length(cartesian_idxs) == 2
#    @assert all(atom_idxs .<= N_atoms) && all(atom_idxs .>= 1) "Atom indexes out of range, must be in 1:$(N_atoms)"
#    @assert all(cartesian_idxs .<= 3) && all(cartesian_idxs .>= 1) "Cartesian indices must be 1, 2, or 3"
   @assert atom_idxs[1] != atom_idxs[2] "Self terms should be calculated with ASR"

#    force_unit = energy_unit(pot) / length_unit(pot)
#    Fᵦ₀ = force_loop_j(pot, posns, force_unit, r_cut, sys_eq.box_sizes_SC, N_atoms, atom_idxs[2], cartesian_idxs[2])
#    posns[atom_idxs[1]][cartesian_idxs[1]] += h
#    ΔFᵦ = force_loop_j(pot, posns, force_unit, r_cut, sys_eq.box_sizes_SC, N_atoms, atom_idxs[2], cartesian_idxs[2])
#    posns[atom_idxs[1]][cartesian_idxs[1]] -= h
#    return - (ΔFᵦ - Fᵦ₀)/h 

    energies = zeros(4)*energy_unit(pot)
    combos = [[h,h],[-h,-h],[h,-h],[-h,h]]

    for (c,combo) in enumerate(combos)

        posns[atom_idxs[1]][cartesian_idxs[1]] += combo[1]
        posns[atom_idxs[2]][cartesian_idxs[2]] += combo[2]

        energies[c] = energy_loop(pot, posns, sys_eq.box_sizes_SC, N_atoms, pot.r_cut)

        #Un-modify
        posns[atom_idxs[1]][cartesian_idxs[1]] -= combo[1]
        posns[atom_idxs[2]][cartesian_idxs[2]] -= combo[2]
    end

    return (1/(4*h^2))*(energies[1] + energies[2] - energies[3] - energies[4])
end