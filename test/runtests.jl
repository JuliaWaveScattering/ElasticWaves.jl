include("../src/ElasticWaves.jl")
using Test

ω = 250000.0;

steel = Elasticity(2; ρ = 7800.0, cp = 5000.0-0.5*im, cs = 3500.0-0.35*im)
bearing = Bearing(2; medium=steel, r1=1.0, r2=2.0)

stress_bcs = [StressBoundary(inner = true), StressBoundary(outer = true)]

M1 = system_matrix(ω, 2, bearing, stress_bcs)

# M2 = stress_matrix_full(2; ω=ω, bearing=bearing)
# M1 - M2 == zeros(Complex{Float64},4,4)

Vector{StressBoundary} <: Vector{AbstractBoundaryCondition}
Vector{StressBoundary} <: (Vector{BC} where BC <: AbstractBoundaryCondition)

# How long does it take to get the coefficients?

basis_order = 121
forcing = rand(basis_order,4) + rand(basis_order,4) .* im
@time mode_coefficients(ω, bearing, stress_bcs, forcing)

# Do both coefficient functions give the same output?
