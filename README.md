# ForceConstants.jl

This package contains code to calculate interatomic force consatnts for pair and three-body potentials. Analytical derivatives up to third order are provided for pair potentials (fastest). Automatic differentiation can also be used for pair potentials, and must be used for three-body potentials. Currently, only 12-6 Lennard-Jones and Stillinger-Weber silicon are implemented and validated. Other potentials could be implemented and used in this framework though.

Due to a bug in the automatic differentiation library, the third order force constants for Stillinger-Weber take a long time to compile. This code is not meant to compete with existing codes like PhonoPy or ALAMODE, rather it accompanies a paper which will be released soon.

The main functions provided by this library are:
- `second_order(sys::SuperCellSystem{D}, pot::Potential, calc::ForceConstantCalculator; n_threads = Threads.n_threads())`
- `third_order(sys::SuperCellSystem{D}, pot::Potential, calc::ForceConstantCalculator; n_threads = Threads.n_threads())`

These return a `DenseForceConstants` object which contains the ifc values, unit and tolerance. These can be accessed with `.values, .units, .tol`. 

-------------------------
There are two kinds of `ForceConstantCalculator`:
- `AnalyticalCalculator(tol, r_cut)`
- `AutoDiffCalculator(tol, r_cut)`
-------------------------
There are two potentials implemented already:
- `LJ(σ, ϵ, r_cut)`
- `StillingerWeberSilicon(; units = true, T = Float64)`
  - The `units` flag controls whether or not Unitful units are used, and `T` controls the precision of the SW parameters
 -------------------------
The last thing required to calculate force constants is a system which contains the atomic positions. The `SuperCellSystem` type has three constructors:
- `SuperCellSystem(positions::AbstractVector{<:AbstractVector}, masses::AbstractVector, box_sizes::AbstractVector, charges = zeros(length(masses))*u"q")`
- `SuperCellSystem(crystal::Crystal{D})`
- `SuperCellSystem(df::DataFrame, masses::AbstractVector, box_sizes::AbstractVector, x_col, y_col, z_col; charges = zeros(length(masses))*u"q")`
-------------------------
### Lennard-Jones Example:
```julia
using ForceConstants
using SimpleCrystals
using Unitful
using StaticArrays

##############
# LJ EXAMPLE #
##############

# Setup system
pot_lj = LJ(3.4u"Å", 0.24037u"kcal * mol^-1", 8.5u"Å")
fcc_crystal = FCC(5.2468u"Å", :Ar, SVector(4,4,4)) #from SimpleCrystals.jl
sys_lj = SuperCellSystem(fcc_crystal);

# Choose force constant calculator
tol = 1e-12
calc_AD_lj = AutoDiffCalculator(tol, pot_lj.r_cut)
calc_analytical_lj = AnalyticalCalculator(tol, pot_lj.r_cut)

# Calculate force constants
ifc2_AD = second_order(sys_lj, pot_lj, calc_AD_lj)
ifc2_analytical = second_order(sys_lj, pot_lj, calc_analytical_lj)
ifc3_AD = third_order(sys_lj, pot_lj, calc_AD_lj)
ifc3_analytical_lj = third_order(sys_lj, pot_lj, calc_analytical_lj)

isapprox(ifc2_analytical.values, ifc2_AD.values, atol = 1e-12)
isapprox(ifc3_analytical.values, ifc3_AD.values, atol = 1e-12)

##############
# SW EXAMPLE #
##############
pot_sw = StillingerWeberSilicon() # energy units are eV / Å^2
diamond_crystal = Diamond(5.43u"Å", :Si, SVector(3,3,3))
sys_sw = SuperCellSystem(diamond_crystal)

tol = 1e-12
calc_AD_sw = AutoDiffCalculator(tol, pot_sw.r_cut)

ifc2_AD_sw = second_order(sys_sw, pot_sw, calc_AD_sw)
ifc3_AD_sw = third_order(sys_sw, pot_sw, calc_AD_sw) # this will take awhile to compile the first time you run it.
```
-------------------------
There are also `second_order!` and `third_order!` which allow you to pass your own storage if you want to re-use memory. A `DenseForceConstants` object will not be returned by these functions, just whatever storage they were provided.
```julia
pot = LJ(3.4u"Å", 0.24037u"kcal * mol^-1", 8.5u"Å")
fcc_crystal = FCC(5.2468u"Å", :Ar, SVector(4,4,4)) #from SimpleCrystals.jl
sys = SuperCellSystem(fcc_crystal);

# Choose force constant calculator
calc_analytical = AnalyticalCalculator(1e-12, pot.r_cut)

# Pre-allocate IFC storage
ifc2 = zeros(Float64, 3*n_atoms(sys), 3*n_atoms(sys))

second_order!(ifc2, sys, pot, calc_analytical) # this modifies `ifc2` in place
```
