@testset "Boundary basis and priors" begin
   

    # Create a basis for the forcing on the inner surface
    ω=4e6
    basis_order = 7;
    numberofsensors = 2
    basis_length = 2*basis_order + 1

    θs = LinRange(0.0, 2pi, basis_length + 1)[1:end-1]
    θ2s = LinRange(0.0, 2pi, 4*basis_length + 1)[1:end-1]
    θs_inv = LinRange(0.0, 2pi, numberofsensors + 1)[1:end-1]
    
    # the pressure and shear fields 
    fp1 = 1*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im
    fs1 = 1*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

    fp2 = 0*exp.(-20.0 .* (θs .- pi/2).^2) + θs .* 0im
    fs2 = 0*exp.(-20.0 .* (θs .- pi/2).^2) + θs .* 0im

    f0 = θs .* 0im

    fouter= 0*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

    # create boundary data for these fields
 
    bd1 = BoundaryData(TractionBoundary(inner=true), θs=θs, fields = hcat(fp1,fouter))
    bd1 = fields_to_fouriermodes(bd1, basis_order)

    bd2 =  BoundaryData(TractionBoundary(outer=true), θs=θs, fields = hcat(fouter, fs1 ))
    bd2 = fields_to_fouriermodes(bd2, basis_order)
    
    bd3 = BoundaryData(TractionBoundary(outer=true), θs=θs, fields = hcat(fp2,fouter))
    bd3 = fields_to_fouriermodes(bd3, basis_order)

    bd4 = BoundaryData(TractionBoundary(outer=true), θs=θs, fields = hcat(fouter,fs2))
    bd4 = fields_to_fouriermodes(bd3, basis_order)


    
    # bd_outer =  BoundaryData(TractionBoundary(inner=true), θs=θs, fields = hcat(θs .* 0im, θs .* 0im))

    boundarybasis1 = BoundaryBasis([bd1,bd2])
    boundarybasis2 = BoundaryBasis([bd3,bd4])



    #FORWARD PROBLEM

    fp=10*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im
    fs=5 *exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

    fouter= 0*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

    bc1_forward = TractionBoundary(inner=true)
    bc2_forward = TractionBoundary(outer=true)

    bc1_inverse = DisplacementBoundary(outer=true)
    bc2_inverse = TractionBoundary(outer=true)


    #bd1_inverse =  BoundaryData(bc1_inverse,θs=θs, fields=hcat(fp,fs))

    bd1_forward =  BoundaryData(bc1_forward,θs=θs, fields=hcat(fp,fs))
    bd1_forward=fields_to_fouriermodes(bd1_forward, basis_order)
    
    bd2_forward = BoundaryData(bc2_forward,θs=θs, fields=hcat(fouter,fouter))
    bd2_forward=fields_to_fouriermodes(bd2_forward, basis_order)


    steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpas = (bearing.outer_radius-bearing.inner_radius ) .* ω ./ steel.cp
    ksas = (bearing.inner_radius-bearing.inner_radius) .* ω ./ steel.cs

    #basis_order = 10
    # basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)
    

    

    #bd1 = BoundaryData(TractionBoundary(inner=true); fourier_modes=forcing_modes[:, 1:2])
    #bd2 = BoundaryData(TractionBoundary(outer=true); fourier_modes=forcing_modes[:, 3:4])

    sim = BearingSimulation(ω, bearing, bd1_forward, bd2_forward)

    wave= ElasticWave(sim)
    #Setting for the inverse problem

    x_outer=[radial_to_cartesian_coordinates([bearing.outer_radius,θ]) for θ in θs_inv ]

    x_inner=[radial_to_cartesian_coordinates([bearing.inner_radius,θ]) for θ in θs]



    traction_outer=[traction(wave,x) for x in x_outer];
    traction_outer=hcat(traction_outer...) |> transpose |>collect



    traction_inner=[traction(wave,x) for x in x_inner]

    traction_inner=hcat(traction_inner...) |> transpose



    #plot(θs_inv,real.(traction_outer[:,1]))
    #plot!(θs_inv,real.(traction_outer[:,2]))


    #plot(θs,real.(traction_inner[:,1]))
    #plot!(θs,real.(traction_inner[:,2]))


    displacement_outer = [displacement(wave,x) for x in x_outer];
    displacement_outer = hcat(displacement_outer...) |> transpose |> collect;

    #plot(θs_inv,real.(displacement_outer[:,1]))
    #plot!(θs_inv,real.(displacement_outer[:,2]))

    bd1_inverse = BoundaryData(
        bc1_inverse;
        θs = θs_inv,
        fields = displacement_outer
    )

    bd1_inverse=fields_to_fouriermodes(bd1_inverse, basis_order)

    bd2_inverse= bd2_forward

 

    inverse_sim = BearingSimulation(ω, bearing, bd1_inverse, bd2_inverse; boundarybasis1=boundarybasis1)
    inverse_sim2 = BearingSimulation(ω, bearing, bd1_inverse, bd2_inverse; boundarybasis1=boundarybasis1, boundarybasis2=boundarybasis2,basis_order=Int(basis_order))
    

   
    inv_wave=ElasticWave(inverse_sim)
    inv_wave2=ElasticWave(inverse_sim2)
   


    errors = maximum(abs.(wave.shear.coefficients - inv_wave.shear.coefficients)) / mean(abs.(wave.shear.coefficients))
    
    errorp=maximum(abs.(wave.pressure.coefficients - inv_wave.pressure.coefficients)) / mean(abs.(wave.pressure.coefficients))

    errors2 = maximum(abs.(wave.shear.coefficients - inv_wave2.shear.coefficients)) / mean(abs.(wave.shear.coefficients))
    
    errorp2=maximum(abs.(wave.pressure.coefficients - inv_wave2.pressure.coefficients)) / mean(abs.(wave.pressure.coefficients))


    
    @test errors < 0.37
    @test errorp < 0.37
    @test errors2 < 0.37
    @test errorp2 < 0.37

    #calculating traction in the inner boundary

    x2_inner = [
        radial_to_cartesian_coordinates([bearing.inner_radius, θ])
    for θ in θ2s]


    x2_outer = [
        radial_to_cartesian_coordinates([bearing.outer_radius, θ])
    for θ in θ2s]


    traction_inv = [traction(inv_wave2,x) for x in x2_inner];
    traction_inv = hcat(traction_inv...) |> transpose |> collect

    traction_inv_out = [traction(inv_wave2,x) for x in x2_outer];
    traction_inv_out = hcat(traction_inv_out...) |> transpose |> collect




    #plot(θs,real.(fp))
    #plot!(θs,real.(fs))

    #plot!(θ2s,real.(traction_inv[:,1]),seriestype=:scatter)
    #plot!(θ2s,real.(traction_inv[:,2]),seriestype=:scatter)

    #savefig("test_prior_method_inner_traction.png")


    #plot(θ2s,real.(traction_inv_out[:,1]),seriestype=:scatter)

    #savefig("test_prior_method_outer_traction.png")


 


    traction_inv2 = [traction(inv_wave,x) for x in x2_inner];
    traction_inv2 = hcat(traction_inv2...) |> transpose |> collect

    traction_inv_out2 = [traction(inv_wave,x) for x in x2_outer];
    traction_inv_out2 = hcat(traction_inv_out2...) |> transpose |> collect




    #plot(θs,real.(fp))
    #plot!(θs,real.(fs))

    #plot!(θ2s,real.(traction_inv2[:,1]),seriestype=:scatter)
    #plot!(θ2s,real.(traction_inv2[:,2]),seriestype=:scatter)

    #savefig("test_prior_method_inner_traction2.png")


    #plot(θ2s,real.(traction_inv_out2[:,1]),seriestype=:scatter)

    #savefig("test_prior_method_outer_traction2.png")


end