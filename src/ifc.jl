export second_order, third_order

function second_order(sys::SuperCellSystem{D}, pot::Potential, 
    calc::ForceConstantCalculator; check_asr = false, check_symmetry = false) where D

    N_atoms = n_atoms(sys)
    Φ = zeros(D*N_atoms,D*N_atoms)

    second_order!(Φ, sys, pot, calc)

    Φ = apply_tols!(Φ, calc.tol)

    check_asr && @assert asr_satisfied(Φ, N_atoms, D, calc.tol) "ASR not satisfied for second order force constants, please report this"
    check_symmetry && @assert issymmetric(Φ) "IFC were not symmetric, please report this"

    return DenseForceConstants(Φ, energy_unit(pot) / length_unit(pot)^2, calc.tol)
end

function third_order(sys::SuperCellSystem{D}, pot::Potential,
    calc::ForceConstantCalculator; check_asr = false) where D

    N_atoms = n_atoms(sys)
    Ψ = zeros(D*N_atoms, D*N_atoms, D*N_atoms)
 
    third_order!(Ψ, sys, pot, calc)

    Ψ = apply_tols!(Ψ, calc.tol)

    Ψ_unit = energy_unit(pot) / length_unit(pot)^3

    check_asr && @assert asr_satisfied(Ψ, N_atoms, D, calc.tol) "ASR not satisfied for third order force constants, please report this"

    return DenseForceConstants(Ψ, Ψ_unit, calc.tol)   
end