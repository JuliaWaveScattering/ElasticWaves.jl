## We define the inverse problem as taking only information from one boundary, such as the outer boundary, and then predicting fields on the other boundary.

@testset "Inverse problems for the modes" begin

    ωs = [1e3,4e4,7e4]

    # NOTE: The inverse problem can become more ill-posed if the imaginary part of the wavespeed is large enogh
    steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0)

    # this non-dimensional number determines what basis_order is neeeded. Also note that kpa[1] and ksa[1] are small, indicating a very low frequency
    kpa = bearing.outer_radius .* ωs ./ steel.cp
    ksa = bearing.outer_radius .* ωs ./ steel.cs

    basis_order = 10
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    # first the forward problem
    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(TractionBoundary(inner = true); fourier_modes = forcing_modes[:,1:2])
    bd2 = BoundaryData(TractionBoundary(outer = true); fourier_modes = forcing_modes[:,3:4])
    
    # Will use only the Fourier modes that were calculated correctly by adding the option: only_stable_modes = true
    method = ModalMethod(tol = 1e-3, only_stable_modes = true)
    sims = [
        BearingSimulation(ω, bearing, bd1, bd2; method = method)
    for ω in ωs];
    waves = ElasticWave.(sims);

    # setup a problem with only boundary information on the outer boundary.
    bd1s = map(waves) do wave
        BoundaryData(TractionBoundary(outer = true);
            fourier_modes = field_modes(wave, bearing.outer_radius, TractionType())
        )
    end

    bd2s = map(waves) do wave
        BoundaryData(DisplacementBoundary(outer = true); 
            fourier_modes = field_modes(wave, bearing.outer_radius, DisplacementType())
        )
    end

    # Calculate the wave only from the outer boundary. We typically think of this as an inverse problem as the (non-zero) traction was specified on the inner boundary
    inverse_sims = [
        BearingSimulation(ωs[i], bearing, bd1s[i], bd2s[i]; method = ModalMethod(only_stable_modes = false))
    for i in eachindex(ωs)];
    inverse_waves = ElasticWave.(inverse_sims);

    # we should recover the first wave
    errors = [
        maximum(abs.(waves[i].potentials[2].coefficients - inverse_waves[i].potentials[2].coefficients)) / mean(abs.(waves[i].potentials[2].coefficients))
    for i in eachindex(ωs)]

    # the lowest frequency here is still unstable
    errors[1] > 1.0
    @test maximum(errors[2:end]) < 1e-13

    errors = [
        maximum(abs.(waves[i].potentials[1].coefficients - inverse_waves[i].potentials[1].coefficients)) / mean(abs.(waves[i].potentials[1].coefficients))
    for i in eachindex(ωs)]

    errors[1] < 1.0
    @test maximum(errors[2:end]) < 1e-13

    # We can redo the lowest frequency with an even more stringent tolerance that excludes unstable modes
    method = ModalMethod(tol = 1e-9, only_stable_modes = true)
    sims = [BearingSimulation(ω, bearing, bd1, bd2; method = method) for ω in ωs];
    waves = ElasticWave.(sims);

    # setup a problem with only boundary information on the outer boundary.
    bd1s = map(waves) do wave
        BoundaryData(TractionBoundary(outer = true);
            fourier_modes = field_modes(wave, bearing.outer_radius, TractionType())
        )
    end

    bd2s = map(waves) do wave
        BoundaryData(DisplacementBoundary(outer = true); 
            fourier_modes = field_modes(wave, bearing.outer_radius, DisplacementType())
        )
    end

    inverse_sims = [
        BearingSimulation(ωs[i], bearing, bd1s[i], bd2s[i]; method = ModalMethod(only_stable_modes = false))
    for i in eachindex(ωs)];
    inverse_waves = ElasticWave.(inverse_sims);

    errors = [
        maximum(abs.(waves[i].potentials[2].coefficients - inverse_waves[i].potentials[2].coefficients)) / mean(abs.(waves[i].potentials[2].coefficients))
    for i in eachindex(ωs)]

    @test maximum(errors) < 1e-13

    errors = [
        maximum(abs.(waves[i].potentials[1].coefficients - inverse_waves[i].potentials[1].coefficients)) / mean(abs.(waves[i].potentials[1].coefficients))
    for i in eachindex(ωs)]

    @test maximum(errors) < 1e-13
    
    # Check if inverse_waves predicts the same traction on the inner boundary
    inner_traction_forcing_modes = [
        field_modes(w, bearing.inner_radius, TractionType())
    for w in inverse_waves]

    # Note how the error in the boundary away from the data is larger
    errors = map(eachindex(ωs)) do i

        in_modes = inner_traction_forcing_modes[i]
        order = Int((size(in_modes)[1] - 1)/2)
        orders = (-order:order) .+ (basis_order + 1)

        maximum(abs.(forcing_modes[orders,1:2] - in_modes)) / mean(abs.(forcing_modes[orders,1:2]))
    end

    @test maximum(errors) < 1e-12

    outer_traction_forcing_modes = [
        field_modes(w, bearing.outer_radius, TractionType())
    for w in inverse_waves]

    errors = map(eachindex(ωs)) do i

        out_modes = outer_traction_forcing_modes[i]
        order = Int((size(out_modes)[1] - 1)/2)
        orders = (-order:order) .+ (basis_order + 1)

        maximum(abs.(forcing_modes[orders,3:4] - out_modes)) / mean(abs.(forcing_modes[orders,3:4]))
    end

    @test maximum(errors) < 1e-13

end

# These tests are incomplete
@testset "Inverse problems for the fields" begin

    ω = 2.0
    steel = Elastic(2; ρ = 7.0, cp = 5.0, cs = 3.5)
    # steel = Elastic(2; ρ = 7.0, cp = 5.0 - 0.1im, cs = 3.5 - 0.06im)
    bearing = RollerBearing(medium = steel, inner_radius = 1.0, outer_radius = 2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpa = bearing.outer_radius * ω / steel.cp
    bearing.inner_radius * ω / steel.cp
    ksa = bearing.outer_radius * ω / steel.cs

    basis_order = 5
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    # 0.0 and 2pi are the same point
    θs = LinRange(0.0,2pi, basis_length + 1)[1:end-1]

    # let's create a focused pressure on the inner boundary
    fs = θs .* 0 + θs .* 0im |> collect;
    fp = θs .* 0 + θs .* 0im |> collect;

    # fp[1:3] = 1e5 .* ones(3) - 1e5im .* ones(3)
    fp[1] = 1e5

    θ0 = θs[1]
    x0 = radial_to_cartesian_coordinates([bearing.inner_radius, θ0])

    bd1 = BoundaryData(TractionBoundary(inner = true); θs = θs, fields = hcat(fp,fs))
    bd2 = BoundaryData(TractionBoundary(outer = true); θs = θs, fields = hcat(fs,fs))

    method = ModalMethod(only_stable_modes = false)
    sim = BearingSimulation(ω, bearing, bd1, bd2; method = method)

    # let's have a look at the modes that were calculated for this BearingSimulation. This is the field we will actual approximate
    inner_field = fouriermodes_to_fields(θs,sim.boundarydata1.fourier_modes)

    @test norm(inner_field[:,1] - fp) < 1e-10

    wave = ElasticWave(sim);
end