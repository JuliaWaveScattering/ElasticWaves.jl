@testset "Boundary basis and priors" begin
   

    # Create a basis for the forcing on the inner surface
    ω=1e6
    basis_order = 3;
    numberofsensors = 2
    basis_length = 2*basis_order + 1

    θs = LinRange(0.0, 2pi, basis_length + 1)[1:end-1]
    θ2s = LinRange(0.0, 2pi, 4*basis_length + 1)[1:end-1]
    θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]
    
    
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

    #ϵ = 0.01*max(bd1_inverse)

    #bd1_inverse= bd1_inverse .+ ϵ.*randn(length(bd1_inverse))

 

    inverse_sim = BearingSimulation(ω, bearing, bd1_inverse, bd2_inverse; boundarybasis1=boundarybasis1)
    inverse_sim2 = BearingSimulation(ω, bearing, bd1_inverse, bd2_inverse; boundarybasis1=boundarybasis1, boundarybasis2=boundarybasis2,basis_order=Int(basis_order))
    

   
    inv_wave=ElasticWave(inverse_sim)
    inv_wave2=ElasticWave(inverse_sim2)
   


    errors = maximum(abs.(wave.potentials[2].coefficients - inv_wave.potentials[2].coefficients)) / mean(abs.(wave.potentials[2].coefficients))
    
    errorp = maximum(abs.(wave.potentials[1].coefficients - inv_wave.potentials[1].coefficients)) / mean(abs.(wave.potentials[1].coefficients))

    errors2 = maximum(abs.(wave.potentials[2].coefficients - inv_wave2.potentials[2].coefficients)) / mean(abs.(wave.potentials[2].coefficients))
    errorp2 = maximum(abs.(wave.potentials[1].coefficients - inv_wave2.potentials[1].coefficients)) / mean(abs.(wave.potentials[1].coefficients))

    
    
    @test errors < 0.37
    @test errorp < 0.37
    @test errors2 < 0.37
    @test errorp2 < 0.37


 


    # the low frequencies are a bit unstable I think due to the hankelh1 singularity
   

    



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


@testset "boudary basis test for the modes with two basis" begin
    
    ωs = [20.0,4e4,7e4]

    # The inverse problem is ill-posed if the wavespeed is complex
    steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0)

    # this non-dimensional number determines what basis_order is neeeded. Also note that kpa[1] and ksa[1] are small, indicating a very low frequency
    kpa = bearing.outer_radius .* ωs ./ steel.cp
    ksa = bearing.outer_radius .* ωs ./ steel.cs

    # estimate the largest basis_order that a wave scattered from the inner boundary can be measured at the outer boundary

    # basis_order = estimate_basisorder(ω,bearing; tol =1e-5)
    basis_order = 12
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    # first the forward problem
    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1_forward = BoundaryData(TractionBoundary(inner = true); fourier_modes = forcing_modes[:,1:2])
    bd2_forward = BoundaryData(TractionBoundary(outer = true); fourier_modes = forcing_modes[:,3:4])

    sims = [BearingSimulation(ω, bearing, bd1_forward, bd2_forward) for ω in ωs];
    waves = ElasticWave.(sims);

    #set priors


    bd1= BoundaryData(TractionBoundary(inner = true); fourier_modes = hcat(forcing_modes[:,1],0*forcing_modes[:,1]))
    bd2= BoundaryData(TractionBoundary(inner = true); fourier_modes = hcat(0*forcing_modes[:,2],forcing_modes[:,2]))

    bbasis1=BoundaryBasis([bd1,bd2])


    bd3= BoundaryData(TractionBoundary(outer  = true); fourier_modes = hcat(forcing_modes[:,3],0*forcing_modes[:,3]))
    bd4= BoundaryData(TractionBoundary(outer  = true); fourier_modes = hcat(0*forcing_modes[:,4],forcing_modes[:,4]))

    bbasis1=BoundaryBasis([bd1,bd2])
    bbasis2=BoundaryBasis([bd3,bd4])

    # setup a problem with only boundary information on the outer boundary
    bd2s = [
        BoundaryData(TractionBoundary(outer = true);
            fourier_modes = field_modes(w, bearing.outer_radius, TractionType())
        )
    for w in waves]

    bd1s = [
        BoundaryData(DisplacementBoundary(outer = true); fourier_modes = field_modes(w, bearing.outer_radius, DisplacementType()))
    for w in waves]

    # calculate the wave from the mixed boundary conditions
    inverse_sims = [
        BearingSimulation(ωs[i], bearing, bd1s[i], bd2s[i],boundarybasis1=bbasis1,boundarybasis2=bbasis2)
    for i in eachindex(ωs)];
    inverse_waves = ElasticWave.(inverse_sims);

    # we should recover the first wave
    errors = [
        maximum(abs.(waves[i].potentials[2].coefficients - inverse_waves[i].potentials[2].coefficients)) / mean(abs.(waves[i].potentials[2].coefficients))
    for i in eachindex(ωs)]

    # the low frequencies works better in prior method
    @test errors[1] < 1e-5
    @test maximum(errors[2:end]) < 1e-4

    errors = [
        maximum(abs.(waves[i].potentials[1].coefficients - inverse_waves[i].potentials[1].coefficients)) / mean(abs.(waves[i].potentials[1].coefficients))
    for i in eachindex(ωs)]

    @test errors[1] < 1e-4
    @test maximum(errors[2:end]) < 1e-4

    # Check if inverse_waves predicts the same traction on the inner boundary
    inner_traction_forcing_modes = [
        field_modes(w, bearing.inner_radius, TractionType())
    for w in inverse_waves]

    # Note how the error in the boundary away from the data is larger
    errors = [
        maximum(abs.(forcing_modes[:,1:2] - inner_traction_forcing_modes[i])) / mean(abs.(forcing_modes[:,1:2]))
    for i in eachindex(ωs)]

    @test errors[1] < 2e-4
    @test maximum(errors[2:end]) < 1e-4

    outer_traction_forcing_modes = [
        field_modes(w, bearing.outer_radius, TractionType())
    for w in inverse_waves]

    errors = [
        maximum(abs.(forcing_modes[:,3:4] - outer_traction_forcing_modes[i])) / mean(abs.(forcing_modes[:,3:4]))
    for i in eachindex(ωs)]

    @test errors[1] < 1e-8
    @test maximum(errors[2:end]) < 1e-14




end

@testset "boudary basis test for the modes with one basis" begin
    
    ωs = [20.0,4e4,7e4]

    # The inverse problem is ill-posed if the wavespeed is complex
    steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0)

    # this non-dimensional number determines what basis_order is neeeded. Also note that kpa[1] and ksa[1] are small, indicating a very low frequency
    kpa = bearing.outer_radius .* ωs ./ steel.cp
    ksa = bearing.outer_radius .* ωs ./ steel.cs

    # estimate the largest basis_order that a wave scattered from the inner boundary can be measured at the outer boundary

    # basis_order = estimate_basisorder(ω,bearing; tol =1e-5)
    basis_order = 12
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    # first the forward problem
    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1_forward = BoundaryData(TractionBoundary(inner = true); fourier_modes = forcing_modes[:,1:2])
    bd2_forward = BoundaryData(TractionBoundary(outer = true); fourier_modes = forcing_modes[:,3:4])

    sims = [BearingSimulation(ω, bearing, bd1_forward, bd2_forward) for ω in ωs];
    waves = ElasticWave.(sims);

    #set priors


    bd1= BoundaryData(TractionBoundary(inner = true); fourier_modes = hcat(forcing_modes[:,1],0*forcing_modes[:,1]))
    bd2= BoundaryData(TractionBoundary(inner = true); fourier_modes = hcat(0*forcing_modes[:,2],forcing_modes[:,2]))

    bbasis1=BoundaryBasis([bd1,bd2])

    # setup a problem with only boundary information on the outer boundary
    bd2s = [
        BoundaryData(TractionBoundary(outer = true);
            fourier_modes = field_modes(w, bearing.outer_radius, TractionType())
        )
    for w in waves]

    bd1s = [
        BoundaryData(DisplacementBoundary(outer = true); fourier_modes = field_modes(w, bearing.outer_radius, DisplacementType()))
    for w in waves]

    # calculate the wave from the mixed boundary conditions
    inverse_sims = [
        BearingSimulation(ωs[i], bearing, bd1s[i], bd2s[i],boundarybasis1=bbasis1)
    for i in eachindex(ωs)];
    inverse_waves = ElasticWave.(inverse_sims);

    # we should recover the first wave
    errors = [
        maximum(abs.(waves[i].potentials[2].coefficients - inverse_waves[i].potentials[2].coefficients)) / mean(abs.(waves[i].potentials[2].coefficients))
    for i in eachindex(ωs)]

    # the low frequencies works better in prior method
    @test errors[1] < 1e-6
    @test maximum(errors[2:end]) < 1e-6

    errors = [
        maximum(abs.(waves[i].potentials[1].coefficients - inverse_waves[i].potentials[1].coefficients)) / mean(abs.(waves[i].potentials[1].coefficients))
    for i in eachindex(ωs)]

    @test errors[1] < 1e-6
    @test maximum(errors[2:end]) < 1e-6

    # Check if inverse_waves predicts the same traction on the inner boundary
    inner_traction_forcing_modes = [
        field_modes(w, bearing.inner_radius, TractionType())
    for w in inverse_waves]

    # Note how the error in the boundary away from the data is larger
    errors = [
        maximum(abs.(forcing_modes[:,1:2] - inner_traction_forcing_modes[i])) / mean(abs.(forcing_modes[:,1:2]))
    for i in eachindex(ωs)]

    @test errors[1] < 2e-4
    @test maximum(errors[2:end]) < 1e-6

    outer_traction_forcing_modes = [
        field_modes(w, bearing.outer_radius, TractionType())
    for w in inverse_waves]

    errors = [
        maximum(abs.(forcing_modes[:,3:4] - outer_traction_forcing_modes[i])) / mean(abs.(forcing_modes[:,3:4]))
    for i in eachindex(ωs)]

    @test errors[1] < 1e-8
    @test maximum(errors[2:end]) < 1e-14




end