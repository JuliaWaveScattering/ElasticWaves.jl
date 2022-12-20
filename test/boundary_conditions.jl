@testset "Boundary conditions" begin

    ω = 4e4

    steel = Elasticity(2; ρ = 7800.0, cp = 5000.0 -0.5im, cs = 3500.0 -0.5im)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius=2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpa = bearing.outer_radius * ω / steel.cp
    ksa = bearing.inner_radius * ω / steel.cs

    basis_order = Int(round(2.0 * max(abs(kpa),abs(ksa)))) + 1
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    @test (Vector{TractionBoundary} <: Vector{AbstractBoundaryCondition}) == false
    @test Vector{TractionBoundary} <: (Vector{BC} where BC <: AbstractBoundaryCondition)

    ## Test the boundary conditions are formulated correctly
    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
    bd1 = BoundaryData(TractionBoundary(inner = true); fourier_modes = forcing_modes[:,1:2])
    bd2 = BoundaryData(TractionBoundary(outer = true); fourier_modes = forcing_modes[:,3:4])

    sim = BearingSimulation(ω, bearing, bd1, bd2)
    wave = ElasticWave(sim);

    θs = -pi:0.1:pi; θs |> length;
    exps = [
        exp(im * θ * m)
    for θ in θs, m = -basis_order:basis_order];

    # the given forcing
    forcing = exps * forcing_modes;
    forcing_inner = forcing[:,1:2];
    forcing_outer = forcing[:,3:4];

    # the predicted forcing
    xs = [bearing.inner_radius .* [cos(θ),sin(θ)] for θ in θs];
    forcing_inner2 = [traction(x,wave) for x in xs];
    forcing_inner2 = hcat(forcing_inner2...) |> transpose;

    xs = [bearing.outer_radius .* [cos(θ),sin(θ)] for θ in θs];
    forcing_outer2 = [traction(x,wave) for x in xs];
    forcing_outer2 = hcat(forcing_outer2...) |> transpose;

    @test maximum(abs.(forcing_inner - forcing_inner2)) < 1e-10
    @test maximum(abs.(forcing_outer - forcing_outer2)) < 1e-10

    ## Test that we recover the displacement given
    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(
        DisplacementBoundary(inner = true);
        fourier_modes = forcing_modes[:,1:2]
    )
    bd2 = BoundaryData(
        DisplacementBoundary(outer = true);
        fourier_modes = forcing_modes[:,3:4]
    )

    sim = BearingSimulation(ω, bearing, bd1, bd2)
    wave = ElasticWave(sim);

    θs = -pi:0.1:pi; θs |> length;
    exps = [exp(im * θ * m) for θ in θs, m = -basis_order:basis_order];

    forcing = exps * forcing_modes;
    forcing_inner = forcing[:,1:2];
    forcing_outer = forcing[:,3:4];

    xs = [bearing.inner_radius .* [cos(θ),sin(θ)] for θ in θs];
    forcing_inner2 = [displacement(x,wave) for x in xs];
    forcing_inner2 = hcat(forcing_inner2...) |> transpose;

    xs = [bearing.outer_radius .* [cos(θ),sin(θ)] for θ in θs];
    forcing_outer2 = [displacement(x,wave) for x in xs];
    forcing_outer2 = hcat(forcing_outer2...) |> transpose;

    @test maximum(abs.(forcing_inner - forcing_inner2)) < 1e-10
    @test maximum(abs.(forcing_outer - forcing_outer2)) < 1e-10

## Test the displacement and tractions equations are correct.

    ## To do this, we use the traction boundary conditions, then predict displacement, and use these predictions as new boundary conditions to predict the traction

    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
    bd1 = BoundaryData(
        TractionBoundary(inner = true);
        fourier_modes = forcing_modes[:,1:2]
    )
    bd2 = BoundaryData(
        TractionBoundary(outer = true);
        fourier_modes = forcing_modes[:,3:4]
    )

    sim = BearingSimulation(ω, bearing, bd1, bd2)
    wave = ElasticWave(sim);

    traction_forcing_modes = hcat(
        traction_modes(bearing.inner_radius, wave),
        traction_modes(bearing.outer_radius, wave)
    )

    # check traction modes are correct
    @test maximum(abs.(traction_forcing_modes - forcing_modes)) / mean(abs.(forcing_modes)) < 1e-12

    # setup a problem with displacement boundary conditions
    bd1 = BoundaryData(
        DisplacementBoundary(inner = true);
        fourier_modes = displacement_modes(bearing.inner_radius, wave)
    )
    bd2 = BoundaryData(
        DisplacementBoundary(outer = true);
        fourier_modes = displacement_modes(bearing.outer_radius, wave)
    )

    sim = BearingSimulation(ω, bearing, bd1, bd2)
    wave2 = ElasticWave(sim);

    # we should exactly recover the first wave
    @test maximum(abs.(wave.shear.coefficients - wave2.shear.coefficients)) /  mean(abs.(wave.shear.coefficients)) < 1e-12

    @test maximum(abs.(wave.pressure.coefficients - wave2.pressure.coefficients)) / mean(abs.(wave.pressure.coefficients)) < 1e-12

    # check that wave2 will predict the same traction on the boundary as wave
    traction_forcing_modes = hcat(
        traction_modes(bearing.inner_radius, wave2),
        traction_modes(bearing.outer_radius, wave2)
    )

    @test maximum(abs.(traction_forcing_modes - forcing_modes)) / mean(abs.(forcing_modes)) < 1e-12

end
