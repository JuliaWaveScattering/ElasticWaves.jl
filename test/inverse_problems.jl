## We define the inverse problem as taking only information from one boundary, such as the outer boundary, and then predicting fields on the other boundary.

@testset "Inverse problems" begin
    ω = 2e4

    steel = Elasticity(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = RollerBearing(medium=steel, r1=1.0, r2=2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpa = bearing.r2 * ω / steel.cp
    bearing.r1 * ω / steel.cp
    ksa = bearing.r2 * ω / steel.cs

    # estimate the largest basis_order that a wave scattered from the inner boundary can be measured at the outer boundary

    basis_order = Int(round(2.0 * max(abs(kpa),abs(ksa)))) + 1
    ms = 0:basis_order
    ratios = [
        abs(hankelh1(m,ksa * bearing.r2) / hankelh1(m,ksa * bearing.r1))
    for m in ms]

    # when the ratio above is too small, we can not resolve that basis_order
    i = findfirst(ratios .< 1e-6)
    basis_order = isnothing(i) ? basis_order : ms[i]
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    # first the forward problem
    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(TractionBoundary(inner = true); modes = forcing_modes[:,1:2])
    bd2 = BoundaryData(TractionBoundary(outer = true); modes = forcing_modes[:,3:4])

    sim = BearingSimulation(ω, bearing, bd1, bd2)
    wave = ElasticWave(sim);

    # setup a problem with only boundary information on the outer boundary
    bd1 = BoundaryData(TractionBoundary(outer = true); modes = traction_modes(bearing.r2, wave))
    bd2 = BoundaryData(DisplacementBoundary(outer = true); modes = displacement_modes(bearing.r2, wave))

    # calculate the wave from the mixed boundary conditions
    sim = BearingSimulation(ω, bearing, bd1, bd2)
    wave2 = ElasticWave(sim);

    # we should recover the first wave
    @test maximum(abs.(wave.shear.coefficients - wave2.shear.coefficients)) / mean(abs.(wave.shear.coefficients)) < 1e-12

    @test maximum(abs.(wave.pressure.coefficients - wave2.pressure.coefficients)) / mean(abs.(wave.pressure.coefficients)) < 1e-12

    # Check if wave2 predicts the same traction on the inner boundary
    inner_traction_forcing_modes = traction_modes(bearing.r1, wave2)

    # Note how the error in the boundary away from the data is larger
    @test maximum(abs.(forcing_modes[:,1:2] - inner_traction_forcing_modes)) / mean(abs.(forcing_modes[:,1:2])) < 2e-6

    outer_traction_forcing_modes = traction_modes(bearing.r2, wave2)

    @test maximum(abs.(forcing_modes[:,3:4] - outer_traction_forcing_modes)) / mean(abs.(forcing_modes[:,3:4])) < 1e-12

end

## Finally, we can use the Jessica Kent method, that treats the TractionBoundary as the forward problem and the DisplacementBoundary as the backward problem. Then we choose what traction modes best match the measured displacement modes. This way we can put restrictions on traction modes, such as specifying contact points on bearings.
