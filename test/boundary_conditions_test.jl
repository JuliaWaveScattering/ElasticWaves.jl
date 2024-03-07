@testset "Boundary conditions" begin

    # choose a low frequency, a mid-frequency, and high frequency
    ωs = [100,4e4,7e4]

    steel = Elastic(2; ρ = 7800.0, cp = 5000.0 -0.5im, cs = 3500.0 -0.5im)
    steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
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

    # could regularise the lowest frequency, though it doesn't appear necessary

    δs = [1e-8,0.0,0.0]
    tol = 1e-3
    sims = map(eachindex(ωs)) do i
        # the option only_stable_modes = false means the method will try to solve for modes which are ill posed 
        method = ModalMethod(regularisation_parameter = δs[i], 
            tol = tol, 
            only_stable_modes = false 
        )
        BearingSimulation(ωs[i], bearing, bd1, bd2; 
            method = method,
            nondimensionalise = true
        ) 
    end

    waves = [ElasticWave(s) for s in sims];

    # check that the predicted modes of the traction 
    # note: could just use TractionType() instead of bd1.boundarytype.fieldtype and bd2.boundarytype.fieldtype
    traction_forcing_modes = [
        hcat(
            field_modes(wave, bearing.inner_radius, bd1.boundarytype.fieldtype),
            field_modes(wave, bearing.outer_radius, bd1.boundarytype.fieldtype)
        )
    for wave in waves]
    
    # it helps to compare the errors mode by mode from n = 0 until n = basis_order 
    traction_errors = [
        sum(abs2,traction_modes - forcing_modes, dims=2)[basis_order:2basis_order+1]
    for traction_modes in traction_forcing_modes]    

    # as the mode increase, or the frequency decreases, the problem becomes more ill-conditioned
    traction_errors[1] |> norm > 1.0

    # however the first modes are reasonably conditioned
    @test traction_errors[1][1:3] |> norm < 8e-3

    # for higher frequencies in this case all the errors are small
    @test traction_errors[2] |> norm < 1e-20
    @test traction_errors[3] |> norm < 1e-20 

    # now we show how to compare the fields in the boundary, instead of the fourier modes
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

    errors[1] > 1e-2
    @test errors[2] < 1e-12
    @test errors[3] < 1e-12

    errors = [maximum(abs.(forcing_outer - f)) for f in forcing_outers];

    errors[1] > 1
    @test errors[2] < 1e-12
    @test errors[3] < 1e-12


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

    # include ill posed modes by setting only_stable_modes = false
    method = ModalMethod(only_stable_modes = false)
    sims = [BearingSimulation(ω, bearing, bd1, bd2; method = method) for ω in ωs];
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

    errors[1] > 0.1
    @test errors[2] < 1e-12
    @test errors[3] < 2e-12

    xs = [bearing.outer_radius .* [cos(θ),sin(θ)] for θ in θs];

    forcing_outers = [
        hcat(([displacement(w, x) for x in xs])...) |> transpose
    for w in waves];

    errors = [maximum(abs.(forcing_outer - f)) for f in forcing_outers];

    errors[1] > 1.0
    @test errors[2] < 1e-10
    @test errors[3] < 1e-10

## Give traction on the boundary then predict the displacement

    ## To do this, we use the traction boundary conditions, then predict displacement, and use these predictions as new boundary conditions to predict the traction

    # need to use a large enough frequency for basis_order = 10
    ω = 1e4

    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
    bd1 = BoundaryData(
        TractionBoundary(inner = true);
        fourier_modes = forcing_modes[:,1:2]
    )
    bd2 = BoundaryData(
        TractionBoundary(outer = true);
        fourier_modes = forcing_modes[:,3:4]
    )

    sim = BearingSimulation(ω, bearing, bd1, bd2; method = method)
    wave = ElasticWave(sim);

    traction_forcing_modes = hcat(
        field_modes(wave, bearing.inner_radius, TractionType()),
        field_modes(wave, bearing.outer_radius, TractionType())
    )

    # check traction modes are correct
    @test maximum(abs.(traction_forcing_modes - forcing_modes)) / mean(abs.(forcing_modes)) < 1e-10

    # setup a problem with displacement boundary conditions
    bd1 = BoundaryData(
        DisplacementBoundary(inner = true);
        fourier_modes = field_modes(wave, bearing.inner_radius, DisplacementType())
    )
    bd2 = BoundaryData(
        DisplacementBoundary(outer = true);
        fourier_modes = field_modes(wave, bearing.outer_radius, DisplacementType())
    )

    sim = BearingSimulation(ω, bearing, bd1, bd2; method = method)
    wave2 = ElasticWave(sim);

    # we should recover the first wave
    @test maximum(abs.(wave.potentials[2].coefficients - wave2.potentials[2].coefficients)) /  mean(abs.(wave.potentials[2].coefficients)) < 1e-10

    @test maximum(abs.(wave.potentials[1].coefficients - wave2.potentials[1].coefficients)) / mean(abs.(wave.potentials[1].coefficients)) < 1e-10

    # check that wave2 will predict the same traction on the boundary as wave
    traction_forcing_modes = hcat(
        field_modes(wave2, bearing.inner_radius, TractionType()),
        field_modes(wave2, bearing.outer_radius, TractionType())
    )

    @test maximum(abs.(traction_forcing_modes - forcing_modes)) / mean(abs.(forcing_modes)) < 1e-10

## Stability check by adding Gaussian noise

    # setup a problem with displacement boundary conditions
    # add 1% error to boundary conditions
    ε = 0.01 * maximum(abs.(field_modes(wave, bearing.inner_radius, DisplacementType())))
    error = ε .* (rand(basis_length, 2) .- 0.5) + ε .* (rand(basis_length, 2) .- 0.5) .* im

    bd1 = BoundaryData(
        DisplacementBoundary(inner=true);
        fourier_modes = field_modes(wave, bearing.inner_radius, DisplacementType()) + error
    )

    ε = 0.01 * maximum(abs.(field_modes(wave, bearing.outer_radius, DisplacementType())))
    error = ε .* (rand(basis_length, 2) .- 0.5) + ε .* (rand(basis_length, 2) .- 0.5) .* im

    bd2 = BoundaryData(
        DisplacementBoundary(outer=true);
        fourier_modes = field_modes(wave, bearing.outer_radius, DisplacementType()) + error
    )

    method = ModalMethod(regularisation_parameter = 0.0, only_stable_modes = false)
    sim = BearingSimulation(ω, bearing, bd1, bd2; method = method)
    wave2 = ElasticWave(sim);

    # the problem is unstable
    maximum(abs.(wave.potentials[2].coefficients - wave2.potentials[2].coefficients)) / mean(abs.(wave.potentials[2].coefficients)) > 0.1
    
    maximum(abs.(wave.potentials[2].coefficients - wave2.potentials[2].coefficients)) / mean(abs.(wave.potentials[2].coefficients)) > 0.1
    
    # the instability is because the problem is ill posed when either decreasing the frequency or increasing the basis_order. Let us keep only_stable_modes to solve this, and only solve for smaller modes

    basis_order = 5
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    # this is still a moderately low frequency
    kpa = (bearing.outer_radius - bearing.inner_radius ) * ω / steel.cp
    ksa = (bearing.outer_radius - bearing.inner_radius) * ω / steel.cs

    forcing_modes = rand(basis_length, 4) + rand(basis_length, 4) .* im
    bd1 = BoundaryData(
        TractionBoundary(inner=true);
        fourier_modes=forcing_modes[:, 1:2]
    )
    bd2 = BoundaryData(
        TractionBoundary(outer=true);
        fourier_modes=forcing_modes[:, 3:4]
    )

    method = ModalMethod(tol = 1e-1, regularisation_parameter = 1e-10, only_stable_modes = true)
    sim = BearingSimulation(ω, bearing, bd1, bd2; method = method)
    wave = ElasticWave(sim);

    # setup a problem with displacement boundary conditions
    # add 1% error to boundary conditions
    ε = 0.01 * maximum(abs.(field_modes(wave, bearing.inner_radius, DisplacementType())));
    error = ε .* (rand(basis_length, 2) .- 0.5) + ε .* (rand(basis_length, 2) .- 0.5) .* im;

    bd1 = BoundaryData(
        DisplacementBoundary(inner=true);
        fourier_modes=field_modes(wave, bearing.inner_radius, DisplacementType()) + error
    );

    ε = 0.01 * maximum(abs.(field_modes(wave, bearing.outer_radius, DisplacementType())))
    error = ε .* (rand(basis_length, 2) .- 0.5) + ε .* (rand(basis_length, 2) .- 0.5) .* im;

    bd2 = BoundaryData(
        DisplacementBoundary(outer=true);
        fourier_modes=field_modes(wave, bearing.outer_radius, DisplacementType()) + error
    );

    sim = BearingSimulation(ω, bearing, bd1, bd2; method = method);
    wave2 = ElasticWave(sim);

    # the problem is still sensitive to errors, but a managable tolerance
    @test maximum(abs.(wave.potentials[2].coefficients - wave2.potentials[2].coefficients)) / mean(abs.(wave.potentials[2].coefficients)) < 0.18
   
    @test maximum(abs.(wave.potentials[1].coefficients - wave2.potentials[1].coefficients)) / mean(abs.(wave.potentials[1].coefficients)) < 0.18

    # but the average error is at most 5 times the noise
    @test mean(abs.(wave.potentials[1].coefficients - wave2.potentials[1].coefficients)) / mean(abs.(wave.potentials[2].coefficients)) < 0.05
    
    @test mean(abs.(wave.potentials[2].coefficients - wave2.potentials[2].coefficients)) / mean(abs.(wave.potentials[2].coefficients)) < 0.05
end
