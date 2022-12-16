include("../src/ElasticWaves.jl")
using Test

ω = 1e5

steel = Elasticity(2; ρ = 7800.0, cp = 5000.0-0.5*im, cs = 3500.0-0.35*im)
bearing = Bearing(2; medium=steel, r1=1.0, r2=2.0)

# this non-dimensional number determines what basis_order is neeeded
kpa = bearing.r2 * ω / steel.cp
ksa = bearing.r1 * ω / steel.cs

# basis_order = 55
basis_order = Int(round(2.0 * max(abs(kpa),abs(ksa)))) + 1
bcs = [TractionBoundary(inner = true), TractionBoundary(outer = true)]

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

@test (Vector{TractionBoundary} <: Vector{AbstractBoundaryCondition}) == false
@test Vector{TractionBoundary} <: (Vector{BC} where BC <: AbstractBoundaryCondition)

# How long does it take to get the coefficients?

# basis_order = 55
basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

## Test that we recover the forcing given

bcs = [TractionBoundary(inner = true), TractionBoundary(outer = true)]

forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
@time wave = ElasticWave(ω, bearing, bcs, forcing_modes);

θs = -pi:0.1:pi; θs |> length;
exps = [exp(im * θ * m) for θ in θs, m = -basis_order:basis_order];

forcing = exps * forcing_modes;
forcing_inner = forcing[:,1:2];
forcing_outer = forcing[:,3:4];

xs = [bearing.r1 .* [cos(θ),sin(θ)] for θ in θs];
forcing_inner2 = [traction(x,wave) for x in xs];
forcing_inner2 = hcat(forcing_inner2...) |> transpose;

xs = [bearing.r2 .* [cos(θ),sin(θ)] for θ in θs];
forcing_outer2 = [traction(x,wave) for x in xs];
forcing_outer2 = hcat(forcing_outer2...) |> transpose;

@test maximum(abs.(forcing_inner - forcing_inner2)) < 1e-10
@test maximum(abs.(forcing_outer - forcing_outer2)) < 1e-10

## Test that we recover the displacement given

bcs = [DisplacementBoundary(inner = true), DisplacementBoundary(outer = true)]

displacement_modes = rand(basis_length,4) + rand(basis_length,4) .* im
@time wave = ElasticWave(ω, bearing, bcs, forcing_modes);

θs = -pi:0.1:pi; θs |> length;
exps = [exp(im * θ * m) for θ in θs, m = -basis_order:basis_order];

forcing = exps * forcing_modes;
forcing_inner = forcing[:,1:2];
forcing_outer = forcing[:,3:4];

xs = [bearing.r1 .* [cos(θ),sin(θ)] for θ in θs];
forcing_inner2 = [displacement(x,wave) for x in xs];
forcing_inner2 = hcat(forcing_inner2...) |> transpose;

xs = [bearing.r2 .* [cos(θ),sin(θ)] for θ in θs];
forcing_outer2 = [displacement(x,wave) for x in xs];
forcing_outer2 = hcat(forcing_outer2...) |> transpose;

@test maximum(abs.(forcing_inner - forcing_inner2)) < 1e-10
@test maximum(abs.(forcing_outer - forcing_outer2)) < 1e-10


## Use the traction boundary conditions, predict displacement, then solve displacement boundary conditions to predict the traction

# bcs = [TractionBoundary(inner = true), TractionBoundary(outer = true)]
#
# forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
#
# # zero the traction on the outer boundary to imitate real bearing
# forcing_modes[:,3:4] = forcing_modes[:,3:4] .* 0
#
# @time wave = ElasticWave(ω, bearing, bcs, forcing_modes);
#
# θs = -pi:0.1:pi; θs |> length;
# exps = [exp(im * θ * m) for θ in θs, m = -basis_order:basis_order];
#
# forcing = exps * forcing_modes;
# forcing_inner = forcing[:,1:2];
# forcing_outer = forcing[:,3:4];
#
# xs = [bearing.r1 .* [cos(θ),sin(θ)] for θ in θs];
# forcing_inner2 = [traction(x,wave) for x in xs];
# forcing_inner2 = hcat(forcing_inner2...) |> transpose;
#
# xs = [bearing.r2 .* [cos(θ),sin(θ)] for θ in θs];
# forcing_outer2 = [traction(x,wave) for x in xs];
# forcing_outer2 = hcat(forcing_outer2...) |> transpose;
#
# bcs = [DisplacementBoundary(inner = true), DisplacementBoundary(outer = true)]
