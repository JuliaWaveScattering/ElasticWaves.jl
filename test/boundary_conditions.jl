@testset "Boundary conditions" begin

    ω = 4e4

    steel = Elasticity(2; ρ = 7800.0, cp = 5000.0 -0.5im, cs = 3500.0 -0.5im)
    bearing = Bearing(2; medium=steel, r1=1.0, r2=2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpa = bearing.r2 * ω / steel.cp
    ksa = bearing.r1 * ω / steel.cs

    basis_order = Int(round(2.0 * max(abs(kpa),abs(ksa)))) + 1
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    bcs = [
        TractionBoundary(inner = true),
        TractionBoundary(outer = true)
    ]

    function big_boundary_system(basis_order)
        map(-basis_order:basis_order) do n
            boundarycondition_system(ω, bearing, bcs, n)
        end
    end

    @time M1s = big_boundary_system(basis_order);

    @test (Vector{TractionBoundary} <: Vector{AbstractBoundaryCondition}) == false
    @test Vector{TractionBoundary} <: (Vector{BC} where BC <: AbstractBoundaryCondition)

    ## Test the boundary conditions are formulated correctly
    bcs = [TractionBoundary(inner = true), TractionBoundary(outer = true)]

    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
    @time wave = ElasticWave(ω, bearing, bcs, forcing_modes);

    θs = -pi:0.1:pi; θs |> length;
    exps = [
        exp(im * θ * m)
    for θ in θs, m = -basis_order:basis_order];

    # the given forcing
    forcing = exps * forcing_modes;
    forcing_inner = forcing[:,1:2];
    forcing_outer = forcing[:,3:4];

    # the predicted forcing
    xs = [bearing.r1 .* [cos(θ),sin(θ)] for θ in θs];
    forcing_inner2 = [traction(x,wave) for x in xs];
    forcing_inner2 = hcat(forcing_inner2...) |> transpose;

    xs = [bearing.r2 .* [cos(θ),sin(θ)] for θ in θs];
    forcing_outer2 = [traction(x,wave) for x in xs];
    forcing_outer2 = hcat(forcing_outer2...) |> transpose;

    @test maximum(abs.(forcing_inner - forcing_inner2)) < 1e-10
    @test maximum(abs.(forcing_outer - forcing_outer2)) < 1e-10

    ## Test that we recover the displacement given
    bcs = [
        DisplacementBoundary(inner = true),
        DisplacementBoundary(outer = true)
    ]

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

## Test the displacement and tractions equations are correct.

    ## To do this, we use the traction boundary conditions, then predict displacement, and use these predictions as new boundary conditions to predict the traction

    bcs = [
        TractionBoundary(inner = true),
        TractionBoundary(outer = true)
    ]

    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
    @time wave = ElasticWave(ω, bearing, bcs, forcing_modes);

    traction_forcing_modes = hcat(
        traction_modes(bearing.r1, wave),
        traction_modes(bearing.r2, wave)
    )

    # check traction modes are correct
    @test maximum(abs.(traction_forcing_modes - forcing_modes)) / maximum(abs.(forcing_modes)) < 1e-10

    # setup a problem with displacement boundary conditions
    displacement_bcs = [
        DisplacementBoundary(inner = true),
        DisplacementBoundary(outer = true)
    ]

    # Calculate the displacement modes on boundaries
    displacement_forcing_modes = hcat(
        displacement_modes(bearing.r1, wave), displacement_modes(bearing.r2, wave)
    )

    # calculate the wave from the displacement boundary conditions
    wave2 = ElasticWave(ω, bearing, displacement_bcs, displacement_forcing_modes);

    # we should exactly recover the first wave
    @test maximum(abs.(wave.shear.coefficients - wave2.shear.coefficients)) /  maximum(abs.(wave.shear.coefficients)) < 1e-11

    @test maximum(abs.(wave.pressure.coefficients - wave2.pressure.coefficients)) / maximum(abs.(wave.pressure.coefficients)) < 1e-11

    # check that wave3 will predict the same traction on the boundary as wave
    traction_forcing_modes = hcat(
        traction_modes(bearing.r1, wave2),
        traction_modes(bearing.r2, wave2)
    )

    @test maximum(abs.(traction_forcing_modes - forcing_modes)) / maximum(abs.(forcing_modes)) < 1e-10

end
