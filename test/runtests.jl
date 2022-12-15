include("../src/ElasticWaves.jl")
using Test

ω = 250000.0;

steel = Elasticity(2; ρ = 7800.0, cp = 5000.0-0.5*im, cs = 3500.0-0.35*im)
bearing = Bearing(2; medium=steel, r1=1.0, r2=2.0)

bcs = [TractionBoundary(inner = true), TractionBoundary(outer = true)]

basis_order = 55

function big_boundary_system(basis_order)
    map(-basis_order:basis_order) do n
        boundarycondition_system(ω, n, bearing, bcs)
    end
end

@time M1s = big_boundary_system(basis_order);
@time Ms = boundarycondition_systems(ω, basis_order, bearing, bcs);

p1_j, p1_h = pressure_traction_modes(ω, bearing.r1, basis_order, steel)
s1_j, s1_h = shear_traction_modes(ω, bearing.r1, basis_order, steel)

p2_j, p2_h = pressure_traction_modes(ω, bearing.r2, basis_order, steel)
s2_j, s2_h = shear_traction_modes(ω, bearing.r2, basis_order, steel)

n = 2basis_order + 1
M2 = vcat(
    hcat(p1_j[:,n],p1_h[:,n],s1_j[:,n],s1_h[:,n]),
    hcat(p2_j[:,n],p2_h[:,n],s2_j[:,n],s2_h[:,n])
);

@test M1s[end] == M2
@test M1s == Ms


# M2 = stress_matrix_full(2; ω=ω, bearing=bearing)
# M1 - M2 == zeros(Complex{Float64},4,4)

@test Vector{TractionBoundary} <: Vector{AbstractBoundaryCondition}
@test Vector{TractionBoundary} <: (Vector{BC} where BC <: AbstractBoundaryCondition)

# How long does it take to get the coefficients?

basis_order = 121
basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)
forcing = rand(basis_length,4) + rand(basis_length,4) .* im
wave = ElasticWave(ω, bearing, bcs, forcing)

# Do both coefficient functions give the same output?
