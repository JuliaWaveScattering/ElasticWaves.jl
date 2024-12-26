using ElasticWaves
using Test, Statistics, LinearAlgebra, MultipleScattering

# need a test to check that the equations for displacement and traction were written correctly. One possible test is to check one version of the integral form of the principal of virtual work. For example, we could check that
# \int_{\mathcal B} \rho_0 \ddot u d V = \int_{\partial \mathcal B} \tau d A
# where \tau is the traction on the surface \mathcal B. So this only works when we specify the traction on the whole surface of \mathcal B.

include("signal_processing_test.jl")

## Tests for waves in free space

include("source_scattering_test.jl")

## Tests for translation Matrices

include("source_test.jl")

## Tests for Cylindrical waves
    # an independent check for the formulas of displacement and traction
    include("traction_displacement_test.jl")

    # tests that the boundary conditions are formed correctly, and uniqueness
    include("boundary_conditions_test.jl")

    # tests that the boundary conditions are formed correctly, and uniqueness
    include("inverse_problems_test.jl")
    include("steel_defect_test.jl")
    
    # tests for the PriorMethod that uses a basis for the solution
    include("boundary_basis_example.jl")
    include("boundary_basis_test.jl")
    
    # tests for the loading profile
    include("loading_profile_tests.jl")