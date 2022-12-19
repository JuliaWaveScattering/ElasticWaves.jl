include("../src/ElasticWaves.jl")

using Test

# need a test to check that the equations for displacement and traction were written correctly


# tests that the boundary conditions are formed correctly, and uniqueness
@time include("boundary_conditions.jl")

# tests that the boundary conditions are formed correctly, and uniqueness
@time include("inverse_problems.jl")
