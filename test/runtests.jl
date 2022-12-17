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
        boundarycondition_system(ω, bearing, bcs, n)
    end
end

@time M1s = big_boundary_system(basis_order);

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

forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
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

bcs = [TractionBoundary(inner = true), TractionBoundary(outer = true)]

forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
@time wave = ElasticWave(ω, bearing, bcs, forcing_modes);

# Calculate the displacement modes on boundaries
displacement_forcing_modes = hcat(
    displacement_modes(bearing.r1, wave), displacement_modes(bearing.r2, wave)
)

# setup a problem with displacement boundary conditions
displacement_bcs = [
    DisplacementBoundary(inner = true),
    DisplacementBoundary(outer = true)
]

# calculate the wave from these displacement boundary conditions
wave2 = ElasticWave(ω, bearing, displacement_bcs, displacement_forcing_modes);

# we should exactly recover the first wave
@test norm(wave.shear.coefficients - wave2.shear.coefficients) / norm(wave.shear.coefficients) < 1e-10

@test norm(wave.pressure.coefficients - wave2.pressure.coefficients) / norm(wave.pressure.coefficients) < 1e-10

# Naturally wave2 will predict the same traction on the boundary as wave
traction_forcing_modes = hcat(
    traction_modes(bearing.r1, wave2), traction_modes(bearing.r2, wave2)
)

@test norm(traction_forcing_modes - forcing_modes) / norm(forcing_modes) < 1e-10

## Finally, we can use the Jessica Kent method, that treats the TractionBoundary as the forward problem and the DisplacementBoundary as the backward problem. Then we choose what traction modes best match the measured displacement modes. This way we can put restrictions on traction modes, such as specifying contact points on bearings.
