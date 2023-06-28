@testset "Boundary conditions" begin

    # choose a low frequency, a mid-frequency, and high frequency
    ωs = [10.0,4e4,7e4]

    steel = Elasticity(2; ρ = 7800.0, cp = 5000.0 -0.5im, cs = 3500.0 -0.5im)
    steel = Elasticity(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpas = bearing.outer_radius .* ωs ./ steel.cp
    ksas = bearing.inner_radius .* ωs ./ steel.cs

    basis_order = 10
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    @test (Vector{BoundaryCondition{TractionType}} <: Vector{BoundaryCondition}) == false
    @test Vector{BoundaryCondition{TractionType}} <: (Vector{BC} where BC <: BoundaryCondition)

# Test the traction boundary conditions are formulated correctly
    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(TractionBoundary(inner=true); fourier_modes=forcing_modes[:, 1:2])
    bd2 = BoundaryData(TractionBoundary(outer=true); fourier_modes=forcing_modes[:, 3:4])

    sims = [BearingSimulation(ω, bearing, bd1, bd2) for ω in ωs];
    waves = [ElasticWave(s) for s in sims];

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
    forcing_inners = [
        hcat(([traction(w, x) for x in xs])...) |> transpose
    for w in waves];

    xs = [bearing.outer_radius .* [cos(θ),sin(θ)] for θ in θs];
    forcing_outers = [
        hcat(([traction(w, x) for x in xs])...) |> transpose
    for w in waves];

    errors = [maximum(abs.(forcing_inner - f)) for f in forcing_inners];

    @test errors[1] < 1e-6
    @test errors[2] < 1e-12
    @test errors[3] < 1e-12

    errors = [maximum(abs.(forcing_outer - f)) for f in forcing_outers];

    @test errors[1] < 1e-6
    @test errors[2] < 1e-12
    @test errors[3] < 1e-12

    ## Stability check by adding Gaussian noise
    
    # add 1% error to boundary conditions
    ε = 0.01 * maximum(abs.(forcing_modes));

    bd1.fourier_modes[:,:] = bd1.fourier_modes[:,:] + ε .* rand(basis_length, 2) + ε .* rand(basis_length, 2) .* im
    bd2.fourier_modes[:,:] = bd2.fourier_modes[:,:] + ε .* rand(basis_length, 2) + ε .* rand(basis_length, 2) .* im

    sims = [BearingSimulation(ω, bearing, bd1, bd2) for ω in ωs]
    waves = [ElasticWave(s) for s in sims]

    errors = map(eachindex(ωs)) do i 

        fs = map(-basis_order:basis_order) do m
            coes = vcat(
                waves[i].pressure.coefficients[:, m+basis_order+1], 
                waves[i].shear.coefficients[:, m+basis_order+1]
            )
            boundarycondition_system(ωs[i], bearing, bd1.boundarytype, bd2.boundarytype, m) * coes
        end
        f = hcat(fs...) |> transpose
        # maximum(abs.(f - hcat(bd1.fourier_modes,bd2.fourier_modes)))
        maximum(abs.(f - forcing_modes))
    end
    
    @test errors[1] < 1.0
    @test errors[2] < 0.04
    @test errors[3] < 0.04

## Test displacement boundary conditions 
    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(
        DisplacementBoundary(inner = true);
        fourier_modes = forcing_modes[:,1:2]
    )
    bd2 = BoundaryData(
        DisplacementBoundary(outer = true);
        fourier_modes = forcing_modes[:,3:4]
    )

    sims = [BearingSimulation(ω, bearing, bd1, bd2) for ω in ωs];
    waves = [ElasticWave(s) for s in sims];

    θs = -pi:0.1:pi; θs |> length;
    exps = [exp(im * θ * m) for θ in θs, m = -basis_order:basis_order];

    forcing = exps * forcing_modes;
    forcing_inner = forcing[:,1:2];
    forcing_outer = forcing[:,3:4];

    xs = [bearing.inner_radius .* [cos(θ),sin(θ)] for θ in θs];

    forcing_inners = [
        hcat(([displacement(w, x) for x in xs])...) |> transpose
    for w in waves];

    errors = [maximum(abs.(forcing_inner - f)) for f in forcing_inners];

    @test errors[1] < 1e-6
    @test errors[2] < 1e-12
    @test errors[3] < 1e-12

    xs = [bearing.outer_radius .* [cos(θ),sin(θ)] for θ in θs];

    forcing_outers = [
        hcat(([displacement(w, x) for x in xs])...) |> transpose
    for w in waves];

    errors = [maximum(abs.(forcing_outer - f)) for f in forcing_outers];

    @test errors[1] < 1e-6
    @test errors[2] < 1e-10
    @test errors[3] < 1e-10

## Give traction on the boundary then predict the displacement

    ## To do this, we use the traction boundary conditions, then predict displacement, and use these predictions as new boundary conditions to predict the traction

    ω = 40.0

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
        field_modes(wave, bearing.inner_radius, TractionType()),
        field_modes(wave, bearing.outer_radius, TractionType())
    )

    # check traction modes are correct

    @test maximum(abs.(traction_forcing_modes - forcing_modes)) / mean(abs.(forcing_modes)) < 1e-8

    # setup a problem with displacement boundary conditions
    bd1 = BoundaryData(
        DisplacementBoundary(inner = true);
        fourier_modes = field_modes(wave, bearing.inner_radius, DisplacementType())
    )
    bd2 = BoundaryData(
        DisplacementBoundary(outer = true);
        fourier_modes = field_modes(wave, bearing.outer_radius, DisplacementType())
    )

    sim = BearingSimulation(ω, bearing, bd1, bd2)
    wave2 = ElasticWave(sim);

    # we should exactly recover the first wave
    @test maximum(abs.(wave.shear.coefficients - wave2.shear.coefficients)) /  mean(abs.(wave.shear.coefficients)) < 1e-8

    @test maximum(abs.(wave.pressure.coefficients - wave2.pressure.coefficients)) / mean(abs.(wave.pressure.coefficients)) < 1e-8

    # check that wave2 will predict the same traction on the boundary as wave
    traction_forcing_modes = hcat(
        field_modes(wave2, bearing.inner_radius, TractionType()),
        field_modes(wave2, bearing.outer_radius, TractionType())
    )

    @test maximum(abs.(traction_forcing_modes - forcing_modes)) / mean(abs.(forcing_modes)) < 1e-8





end
