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
    modes = -basis_order:basis_order
    # basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)
    basis_length = length(modes)

    @test (Vector{BoundaryCondition{TractionType}} <: Vector{BoundaryCondition}) == false
    @test Vector{BoundaryCondition{TractionType}} <: (Vector{BC} where BC <: BoundaryCondition)

# Test the traction boundary conditions are formulated correctly
    forcing_coefficients = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(TractionBoundary(inner=true); coefficients = forcing_coefficients[:, 1:2], modes = modes)
    bd2 = BoundaryData(TractionBoundary(outer=true); coefficients = forcing_coefficients[:, 3:4], modes = modes)

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
    traction_mode_errors = [
        sum.(abs2,field_modes(wave, bearing.inner_radius, bd1.boundarytype.fieldtype) - bd1.coefficients) + 
        sum.(abs2,field_modes(wave, bearing.outer_radius, bd2.boundarytype.fieldtype) - bd2.coefficients)
    for wave in waves];

    # as the mode increase, or the frequency decreases, the problem becomes more ill-conditioned
    traction_mode_errors[1] |> norm > 1.0

    # however the first modes are reasonably conditioned
    @test traction_mode_errors[1][1:3,:] |> norm < 8e-3

    # for higher frequencies in this case all the errors are small
    @test traction_mode_errors[2] |> norm < 1e-20
    @test traction_mode_errors[3] |> norm < 1e-20 

    # now we show how to compare the fields in the boundary, instead of the fourier modes
    θs = -pi:0.1:pi; θs |> length;
    exps = [exp(im * θ * m) for θ in θs, m = modes];

    # the given forcing
    forcing = exps * forcing_coefficients;
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
    forcing_coefficients = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(
        DisplacementBoundary(inner = true);
        coefficients = forcing_coefficients[:,1:2],
        modes = modes
    )
    bd2 = BoundaryData(
        DisplacementBoundary(outer = true);
        coefficients = forcing_coefficients[:,3:4],
        modes = modes
    )

    # include ill posed modes by setting only_stable_modes = false
    method = ModalMethod(only_stable_modes = false)
    sims = [BearingSimulation(ω, bearing, bd1, bd2; method = method) for ω in ωs];
    waves = [ElasticWave(s) for s in sims];

    exps = [exp(im * θ * m) for θ in θs, m = modes];

    forcing = exps * forcing_coefficients;
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

    forcing_coefficients = rand(basis_length,4) + rand(basis_length,4) .* im
    bd1 = BoundaryData(
        TractionBoundary(inner = true);
        coefficients = forcing_coefficients[:,1:2],
        modes = modes
    )
    bd2 = BoundaryData(
        TractionBoundary(outer = true);
        coefficients = forcing_coefficients[:,3:4],
        modes = modes
    )

    sim = BearingSimulation(ω, bearing, bd1, bd2; method = method)
    wave = ElasticWave(sim);

    traction_modes_error = norm(field_modes(wave, bearing.inner_radius, TractionType()) - bd1.coefficients) 
    traction_modes_error = traction_modes_error + norm(field_modes(wave, bearing.outer_radius, TractionType()) - bd2.coefficients)

    # check traction modes are correct
    @test traction_modes_error / mean(abs.(forcing_coefficients)) < 1e-10

    # setup a problem with displacement boundary conditions
    bd1_inverse = BoundaryData(
        DisplacementBoundary(inner = true);
        coefficients = field_modes(wave, bearing.inner_radius, DisplacementType()),
        modes = wave.method.modes # need to give correct order of modes! So using modes = modes would be wrong here
    )
    bd2_inverse = BoundaryData(
        DisplacementBoundary(outer = true);
        coefficients = field_modes(wave, bearing.outer_radius, DisplacementType()),
        modes = wave.method.modes
    )

    sim = BearingSimulation(ω, bearing, bd1, bd2; method = method)
    wave_inverse = ElasticWave(sim);

    # we should recover the first wave
    @test maximum(abs.(wave.potentials[2].coefficients - wave_inverse.potentials[2].coefficients)) /  mean(abs.(wave.potentials[2].coefficients)) < 1e-10

    @test maximum(abs.(wave.potentials[1].coefficients - wave_inverse.potentials[1].coefficients)) / mean(abs.(wave.potentials[1].coefficients)) < 1e-10


    # check that wave_inverse will predict the same traction on the boundary as wave
    traction_modes_error = norm(field_modes(wave_inverse, bearing.inner_radius,  bd1.boundarytype.fieldtype) - bd1.coefficients) 
    traction_modes_error = traction_modes_error + norm(field_modes(wave_inverse, bearing.outer_radius, bd2.boundarytype.fieldtype) - bd2.coefficients)

    # check traction modes are correct
    @test traction_modes_error / mean(abs.(forcing_coefficients)) < 1e-10

## Stability check by adding Gaussian noise

    # setup a problem with displacement boundary conditions
    # add 1% error to boundary conditions
    ε = 0.01 * maximum(abs.(field_modes(wave, bearing.inner_radius, DisplacementType())))
    error = ε .* (rand(basis_length, 2) .- 0.5) + ε .* (rand(basis_length, 2) .- 0.5) .* im

    bd1_inverse = BoundaryData(
        DisplacementBoundary(inner=true);
        coefficients = field_modes(wave, bearing.inner_radius, DisplacementType()) + error,
        modes = wave.method.modes
    )

    ε = 0.01 * maximum(abs.(field_modes(wave, bearing.outer_radius, DisplacementType())))
    error = ε .* (rand(basis_length, 2) .- 0.5) + ε .* (rand(basis_length, 2) .- 0.5) .* im

    bd2_inverse = BoundaryData(
        DisplacementBoundary(outer=true);
        coefficients = field_modes(wave, bearing.outer_radius, DisplacementType()) + error,
        modes = wave.method.modes
    )
    
    method = ModalMethod(regularisation_parameter = 0.0, only_stable_modes = false)
    sim = BearingSimulation(ω, bearing, bd1_inverse, bd2_inverse; method = method);
    wave_inverse = ElasticWave(sim);

    # the problem is unstable
    @test maximum(abs.(wave.potentials[2].coefficients - wave_inverse.potentials[2].coefficients)) /  mean(abs.(wave.potentials[2].coefficients)) > 0.1

    @test maximum(abs.(wave.potentials[1].coefficients - wave_inverse.potentials[1].coefficients)) / mean(abs.(wave.potentials[1].coefficients)) > 0.1
  
    # the instability is because the problem is ill posed when either decreasing the frequency or increasing the basis_order. Let us keep only_stable_modes to solve this, and only solve for smaller modes

    basis_order = 5
    modes = -basis_order:basis_order;
    basis_length = length(modes)

    # this is still a moderately low frequency
    kpa = (bearing.outer_radius - bearing.inner_radius ) * ω / steel.cp
    ksa = (bearing.outer_radius - bearing.inner_radius) * ω / steel.cs

    forcing_coefficients = rand(basis_length, 4) + rand(basis_length, 4) .* im
    bd1 = BoundaryData(
        TractionBoundary(inner=true);
        coefficients = forcing_coefficients[:, 1:2],
        modes = modes
    )
    bd2 = BoundaryData(
        TractionBoundary(outer=true);
        coefficients = forcing_coefficients[:, 3:4],
        modes = modes
    )

    method = ModalMethod(tol = 1e-1, regularisation_parameter = 1e-10, only_stable_modes = true)
    sim = BearingSimulation(ω, bearing, bd1, bd2; method = method)
    wave = ElasticWave(sim);

    # setup a problem with displacement boundary conditions
    # add 1% error to boundary conditions
    ε = 0.01 * maximum(abs.(field_modes(wave, bearing.inner_radius, DisplacementType())));
    error = ε .* (rand(basis_length, 2) .- 0.5) + ε .* (rand(basis_length, 2) .- 0.5) .* im;

    bd1_inverse = BoundaryData(
        DisplacementBoundary(inner=true);
        coefficients=field_modes(wave, bearing.inner_radius, DisplacementType()) + error,
        modes = wave.method.modes
    );

    ε = 0.01 * maximum(abs.(field_modes(wave, bearing.outer_radius, DisplacementType())))
    error = ε .* (rand(basis_length, 2) .- 0.5) + ε .* (rand(basis_length, 2) .- 0.5) .* im;

    bd2_inverse = BoundaryData(
        DisplacementBoundary(outer=true);
        coefficients = field_modes(wave, bearing.outer_radius, DisplacementType()) + error,
        modes = wave.method.modes
    );

    sim = BearingSimulation(ω, bearing, bd1_inverse, bd2_inverse; method = method);
    wave_inverse = ElasticWave(sim);

    # the problem is still sensitive to errors, but a managable tolerance
    @test maximum(abs.(wave.potentials[2].coefficients - wave_inverse.potentials[2].coefficients)) / mean(abs.(wave.potentials[2].coefficients)) < 0.18
   
    @test maximum(abs.(wave.potentials[1].coefficients - wave_inverse.potentials[1].coefficients)) / mean(abs.(wave.potentials[1].coefficients)) < 0.18

    # but the average error is at most 5 times the noise
    @test mean(abs.(wave.potentials[1].coefficients - wave_inverse.potentials[1].coefficients)) / mean(abs.(wave.potentials[2].coefficients)) < 0.05
    
    @test mean(abs.(wave.potentials[2].coefficients - wave_inverse.potentials[2].coefficients)) / mean(abs.(wave.potentials[2].coefficients)) < 0.05
end

# test using non centred modes
@testset "Boundary conditions any modes" begin

    # choose a mid-frequency and high frequency
    ωs = [20.0,50.0]

    steel = Elastic(2; ρ = 7.0, cp = 5.0, cs = 3.0)
    bearing = RollerBearing(medium=steel, inner_radius=1.5, outer_radius = 2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpas = bearing.outer_radius .* ωs ./ steel.cp
    ksas = bearing.inner_radius .* ωs ./ steel.cs

    @test (Vector{BoundaryCondition{TractionType}} <: Vector{BoundaryCondition}) == false
    @test Vector{BoundaryCondition{TractionType}} <: (Vector{BC} where BC <: BoundaryCondition)

# Test the traction boundary conditions are formulated correctly
    modes = [-8,-5,1,2,3]
    basis_length = length(modes)
    forcing_coefficients = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(TractionBoundary(inner=true); 
        coefficients = forcing_coefficients[:, 1:2],
        modes = modes
    )
    bd2 = BoundaryData(TractionBoundary(outer=true); 
        coefficients = forcing_coefficients[:, 3:4],
        modes = modes
    )

    tol = 1e-5
    method = ModalMethod(tol = tol)

    sims = map(ωs) do ω
        BearingSimulation(ω, bearing, bd1, bd2; 
            method = method
        ) 
    end

    waves = [ElasticWave(s) for s in sims];

    # check that the predicted modes of the traction 
    # note: could just use TractionType() instead of bd1.boundarytype.fieldtype and bd2.boundarytype.fieldtype
    traction_modes_error = [
        norm(field_modes(wave, bearing.inner_radius,  bd1.boundarytype.fieldtype) - bd1.coefficients) + 
        norm(field_modes(wave, bearing.outer_radius, bd2.boundarytype.fieldtype) - bd2.coefficients)
    for wave in waves]

    # for higher frequencies in this case all the errors are small
    @test traction_modes_error[1] < 1e-12
    @test traction_modes_error[2] < 1e-12

    # Let's repeat, but give some modes to the boundary which can not be solved

    modes = [-8,-18,1,35,3,4,5,6] # -18 is too high for ωs[1] and 35 is too high for ωs[2]
    basis_length = length(modes)
    forcing_coefficients = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(TractionBoundary(inner=true); 
        coefficients = forcing_coefficients[:, 1:2],
        modes = modes
    )
    bd2 = BoundaryData(TractionBoundary(outer=true); 
        coefficients = forcing_coefficients[:, 3:4],
        modes = modes
    )

    tol = 1e-5
    method = ModalMethod(tol = tol)

    sims = map(ωs) do ω
        BearingSimulation(ω, bearing, bd1, bd2; method = method) 
    end;

    waves = [ElasticWave(s) for s in sims];

    # now as not all modes were recovered, we need to select the correct modes from the original data to compare with
    traction_errors = map(waves) do wave

        inds = [findfirst(bd1.modes .== m) for m in wave.method.modes]
        err = sum(abs2,field_modes(wave, bearing.inner_radius, bd1.boundarytype.fieldtype) - bd1.coefficients[inds,:], dims=2)
        
        inds = [findfirst(bd2.modes .== m) for m in wave.method.modes]
        err = err + sum(abs2,field_modes(wave, bearing.outer_radius, bd2.boundarytype.fieldtype) - bd2.coefficients[inds,:], dims=2)

        err
    end

    # for higher frequencies in this case all the errors are small
    @test traction_errors[1] |> maximum < 1e-20
    @test traction_errors[2] |> maximum < 1e-20 
end
