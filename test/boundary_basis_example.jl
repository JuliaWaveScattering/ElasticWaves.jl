# This is very commented example of how to use the PriorMethod

# Here is an example of using a traction basis with three elements that should give an exact inversion
@testset "Boundary basis example" begin

    # the higher the frequency, the worse the result. This is already a high frequency.
    ω = 1e6
    steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0)
    
    dr = bearing.outer_radius - bearing.inner_radius
    kp = (ω / steel.cp)
    kp * dr

    μ = steel.ρ * steel.cs^2
    λ2μ = steel.cp^2 * steel.ρ

## Forward problem with forcing on inner boundary to create a basis
    basis_order = 30;    
    θs = LinRange(0.0, 2pi, 2basis_order+2)[1:end-1]
    # using 2basis_order + 2 guarantees that we can exactly represent the above with Fourier modes with basis_order number of modes

    # choose a basis for the pressure and shear on the inner boundary
    # with twp boundaries we need only one sensor position that measures both pressure and shear
    # θos = [0.0,pi/3,pi];
    θos = [0.0,pi/3];
    fp1s = [
        μ .* exp.(-20.0 .* (sin.(θs ./ 2.0) .- sin(θ ./ 2.0)).^2) + θs .* 0im 
    for θ = θos]; 
        
    fs1s = [
        μ .* 0.5 .* exp.(-20.0 .* (sin.(θs ./ 2.0) .- sin(θ ./ 2.0)).^2) + θs .* 0im 
    for θ = θos];

## The forward problem

    # choose one combination of the basis to be the true case.
    fp1 = sum(fp1s);
    fs1 = sum(fs1s);

    # choose just zero traction for the outer boundary
    fp2 = 0.0 .* fp1 + 0im .* fp1;
    fs2 = 0.0 .* fp1 + 0im .* fp1;

    # Create boundary data for the forward problem
        bc1_forward = TractionBoundary(inner=true)
        bc2_forward = TractionBoundary(outer=true)

        bd1_for = BoundaryData(bc1_forward, θs=θs, fields = hcat(fp1,fs1))
        bd1_for = fields_to_fouriermodes(bd1_for)
    
        # a quick test that the fields and Fourier modes are exactly invertable for convenience 
        bd_test = fouriermodes_to_fields(bd1_for)
        @test norm(bd_test.fields - bd1_for.fields) / norm(bd1_for.fields) < 1e-14
        
        bd2_for = BoundaryData(bc2_forward, θs=θs, fields = hcat(fp2,fs2))

    # Solve the whole field for the forward problem        
        # the method specifies to use only stable modes.
        modal_method = ModalMethod(tol = 1e-6, only_stable_modes = true)
        forward_sim = BearingSimulation(ω, bearing, bd1_for, bd2_for; 
            method = modal_method,
            nondimensionalise = true);
        wave = ElasticWave(forward_sim);

    ## The inverse problem
        bc1_inv = DisplacementBoundary(outer=true)
        bc2_inv = TractionBoundary(outer=true)

        # The traction and displacement field on the outside boundary are similar, because traction for steel is far higher than displacement
        
        # bd1_outer = BoundaryData(bc1_inv, bearing.outer_radius, θs, wave)
        # bd2_outer = BoundaryData(bc2_inv, bearing.outer_radius, θs, wave)
        # using Plots
        # plot(θs, abs.(bd2_outer.fields))
        # plot!(θs, abs.(bd1_outer.fields), linestyle= :dash)

    
    # as there are two basis functions, we will need to at least two measurements. Typically each sensor, or point on the boundary, gives two measurements. So one sensor can be enough. 
        numberofsensors = Int(ceil(length(θos)/2))
        
        θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]

        # create the data from evaluating the forward problem 
        bd1_inverse = BoundaryData(bc1_inv, bearing.outer_radius, θs_inv, wave)
        
    # the traction field is known to be zero everywhere. There needs to be enough points sampled here to match the basis_order used, I think. So instead just using the modal representation that was used for the forward problem.

        bd2_inverse = bd2_for # this is a bit of an inverse crime.

    # Create a boundary basis for the inverse problem for the inner boundary
        bd1s = [
            BoundaryData(bc1_forward, θs = θs, fields = hcat(fp1s[j],fs1s[j]))
        for j in eachindex(fp1s)]

        boundarybasis1 = BoundaryBasis(bd1s)

    # solve the inverse problem with the PriorMethod
    method = PriorMethod(tol = modal_method.tol, modal_method = modal_method)

    inverse_sim = BearingSimulation(ω, method, bearing, bd1_inverse, bd2_inverse;
        boundarybasis1 = boundarybasis1,
        nondimensionalise = true
    );

    inverse_wave = ElasticWave(inverse_sim);

## Test how well we recover the inner traction
    # using a convenience function
    bd1_inner, bd2_outer = boundary_data(forward_sim, inverse_wave)

    # even excluding some modes the recovery is very good with just one sensor
    # this is because the higher order modes (in this case) contribute very little to the field
    @test norm(bd1_inner.fields - bd1_for.fields) / norm(bd1_for.fields) < 1e-10
    
    @test norm(wave.potentials[1].coefficients - inverse_wave.potentials[1].coefficients) / norm(wave.potentials[1].coefficients) < 2e-12

    @test norm(wave.potentials[2].coefficients - inverse_wave.potentials[2].coefficients) / norm(wave.potentials[2].coefficients) < 2e-12

    modes1 = field_modes(wave, bearing.inner_radius, bc1_forward.fieldtype);
    inverse_modes1 = field_modes(inverse_wave, bearing.inner_radius, bc1_forward.fieldtype);

    @test norm(modes1 - inverse_modes1) / norm(modes1) < 2e-12
end
