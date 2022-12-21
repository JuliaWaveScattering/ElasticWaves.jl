## We define the inverse problem as taking only information from one boundary, such as the outer boundary, and then predicting fields on the other boundary.

@testset "Inverse problems for the modes" begin
    ω = 2e4

    steel = Elasticity(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius=2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpa = bearing.outer_radius * ω / steel.cp
    bearing.inner_radius * ω / steel.cp
    ksa = bearing.outer_radius * ω / steel.cs

    # estimate the largest basis_order that a wave scattered from the inner boundary can be measured at the outer boundary

    basis_order = estimate_basisorder(ω,bearing; tol =1e-5)
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    # first the forward problem
    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(TractionBoundary(inner = true); fourier_modes = forcing_modes[:,1:2])
    bd2 = BoundaryData(TractionBoundary(outer = true); fourier_modes = forcing_modes[:,3:4])

    sim = BearingSimulation(ω, bearing, bd1, bd2)
    wave = ElasticWave(sim);

    # setup a problem with only boundary information on the outer boundary
    bd1 = BoundaryData(TractionBoundary(outer = true); fourier_modes = traction_modes(bearing.outer_radius, wave))
    bd2 = BoundaryData(DisplacementBoundary(outer = true); fourier_modes = displacement_modes(bearing.outer_radius, wave))

    # calculate the wave from the mixed boundary conditions
    sim = BearingSimulation(ω, bearing, bd1, bd2)
    wave2 = ElasticWave(sim);

    # we should recover the first wave
    @test maximum(abs.(wave.shear.coefficients - wave2.shear.coefficients)) / mean(abs.(wave.shear.coefficients)) < 1e-12

    @test maximum(abs.(wave.pressure.coefficients - wave2.pressure.coefficients)) / mean(abs.(wave.pressure.coefficients)) < 1e-12

    # Check if wave2 predicts the same traction on the inner boundary
    inner_traction_forcing_modes = traction_modes(bearing.inner_radius, wave2)

    # Note how the error in the boundary away from the data is larger
    @test maximum(abs.(forcing_modes[:,1:2] - inner_traction_forcing_modes)) / mean(abs.(forcing_modes[:,1:2])) < 2e-6

    outer_traction_forcing_modes = traction_modes(bearing.outer_radius, wave2)

    @test maximum(abs.(forcing_modes[:,3:4] - outer_traction_forcing_modes)) / mean(abs.(forcing_modes[:,3:4])) < 1e-12

end

@testset "Inverse problems for the fields" begin
    ω = 2.0
    steel = Elasticity(2; ρ = 7.0, cp = 5.0 - 0.1im, cs = 3.5 - 0.06im)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius=2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpa = bearing.outer_radius * ω / steel.cp
    bearing.inner_radius * ω / steel.cp
    ksa = bearing.outer_radius * ω / steel.cs

    # estimate the largest basis_order that a wave scattered from the inner boundary can be measured at the outer boundary

    basis_order = estimate_basisorder(ω,bearing; tol =1e-5)
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

    sim = BearingSimulation(ω, bearing, bd1, bd2)

    # let's have a look at the modes that were calculated during the Bearing. This is the field we will actual approximate
    inner_field = fouriermodes_to_fields(θs,sim.boundarydata1.fourier_modes)

    @test norm(inner_field[:,1] - fp) < 1e-10

    wave = ElasticWave(sim);
end


## Finally, we can use the Jessica Kent method, that treats the TractionBoundary as the forward problem and the DisplacementBoundary as the backward problem. Then we choose what traction modes best match the measured displacement modes. This way we can put restrictions on traction modes, such as specifying contact points on bearings.
