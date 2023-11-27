@testset "Boundary basis and priors" begin
    
    # Create a basis for the forcing on the inner surface

    basis_order = 3;
    numberofsensors = 10
    basis_length = 2*basis_order + 1

    θs = LinRange(0.0, 2pi, basis_length + 1)[1:end-1]
    
    # the pressure and shear fields 
    fp1 = exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im
    fs1 = exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

    fp2 = exp.(-20.0 .* (θs .- pi/2).^2) + θs .* 0im
    fs2 = exp.(-20.0 .* (θs .- pi/2).^2) + θs .* 0im

    f0 = θs .* 0im

    # create boundary data for these fields

    bd1 = BoundaryData(TractionBoundary(inner=true), θs=θs, fields = hcat(fp1,fs1))
    bd1 = fields_to_fouriermodes(bd1, basis_order)

    bd2 =  BoundaryData(TractionBoundary(inner=true), θs=θs, fields = hcat(fp2,fs2))
    bd2 = fields_to_fouriermodes(bd2, basis_order)
    
    # bd_outer =  BoundaryData(TractionBoundary(inner=true), θs=θs, fields = hcat(θs .* 0im, θs .* 0im))

    boundarybasis1 = BoundaryBasis([bd1,bd2])
    boundarybasis2 = BoundaryBasis([bd1,bd2])



#FORWARD PROBLEM

fp=10*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im
fs=5 *exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

fouter= 0*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im


bd1_inverse =  BoundaryData(bc1_inverse,θs=θs, fields=hcat(fp,fs))

bd1_forward =  BoundaryData(bc1_forward,θs=θs, fields=hcat(fp,fs))
bd2_forward = BoundaryData(bc2_forward,θs=θs, fields=hcat(fouter,fouter))

   
   
    ω = 4e4

    steel = Elasticity(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpas = bearing.outer_radius .* ω ./ steel.cp
    ksas = bearing.inner_radius .* ω ./ steel.cs

    basis_order = 10
    # basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)
    basis_length = 2*basis_order + 1

    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

    bd1 = BoundaryData(TractionBoundary(inner=true); fourier_modes=forcing_modes[:, 1:2])
    bd2 = BoundaryData(TractionBoundary(outer=true); fourier_modes=forcing_modes[:, 3:4])

    sims = [BearingSimulation(ω, bearing, bd1, bd2)]
    wave = ElasticWave.(sims)

end