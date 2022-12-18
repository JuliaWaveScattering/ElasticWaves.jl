## We define the inverse problem as taking only information from one boundary, such as the outer boundary, and then predicting fields on the other boundary.

@testset "Inverse problems" begin
    ω = 2e4

    steel = Elasticity(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = Bearing(2; medium=steel, r1=1.0, r2=2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpa = bearing.r2 * ω / steel.cp
    ksa = bearing.r1 * ω / steel.cs

    # basis_order = 55
    basis_order = Int(round(2.0 * max(abs(kpa),abs(ksa)))) + 1

    # first the forward problem
    bcs = [
        TractionBoundary(inner = true),
        TractionBoundary(outer = true)
    ]

    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
    @time wave = ElasticWave(ω, bearing, bcs, forcing_modes);

    # setup a problem with only boundary information on the outer boundary
    mixed_bcs = [
        TractionBoundary(outer = true),
        DisplacementBoundary(outer = true)
    ]

    # Calculate the boundary data
    mixed_forcing_modes = hcat(
        traction_modes(bearing.r2, wave), displacement_modes(bearing.r2, wave)
    )

    # calculate the wave from the mixed boundary conditions
    wave2 = ElasticWave(ω, bearing, mixed_bcs, mixed_forcing_modes);

    # we should exactly recover the first wave
    @test norm(wave.shear.coefficients - wave2.shear.coefficients) / norm(wave.shear.coefficients) < 1e-10

    maximum(abs.(wave.shear.coefficients - wave2.shear.coefficients)) / maximum(abs.(wave.shear.coefficients))

    @test norm(wave.pressure.coefficients - wave2.pressure.coefficients) / norm(wave.pressure.coefficients) < 1e-10

    maximum(abs.(wave.pressure.coefficients - wave2.pressure.coefficients)) / maximum(abs.(wave.pressure.coefficients))

    # But wave2 does not predict the same traction on the boundary as wave
    traction_forcing_modes2 = hcat(
        traction_modes(bearing.r1, wave2),
        traction_modes(bearing.r2, wave2)
    )

    norm(traction_forcing_modes - traction_forcing_modes2) / norm(traction_forcing_modes)

end

## Finally, we can use the Jessica Kent method, that treats the TractionBoundary as the forward problem and the DisplacementBoundary as the backward problem. Then we choose what traction modes best match the measured displacement modes. This way we can put restrictions on traction modes, such as specifying contact points on bearings.
