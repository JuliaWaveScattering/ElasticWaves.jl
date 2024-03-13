# This is very commented example of how to use the PriorMethod

@testset "Loading profile" begin

# the higher the frequency, the worse the result. This is already a high frequency.
medium = Elastic(2; ρ = 2.0, cp = 1.0 - 0.0im, cs = 0.8 - 0.0im)

Ω = 2pi * 15 / 60 # large wind turbines rotate at about 15 rpm

bearing = RollerBearing(medium = medium, 
    inner_radius = 1.5, outer_radius = 2.0, 
    angular_speed = Ω,  
    rollers_inside = true
)

frequency_order = 3

ωms = natural_frequencies(bearing, frequency_order) |> collect

ω = ωms[end]

dr = bearing.outer_radius - bearing.inner_radius
kp = (ω / medium.cp)
kp * dr

# create the true loading profile, then solve the forward problem to create dat for tthe inverse problem

    loading_resolution = 50;
    loading_θs = LinRange(0.0, 2pi, 2*loading_resolution+2)[1:end-1]

    θo = 3pi/2;
    # fp_loading = 0.2 .- exp.(-0.4 .* (sin.(loading_θs) .- sin(θo)).^2) + loading_θs .* 0im; 
    # fs_loading = 0.1 .* fp_loading;
    fp_loading = 0.0 .+ 0.3 .* cos.(loading_θs) .- 0.1 .* cos.(2 .* loading_θs);
    fs_loading = 0.0 .* fp_loading;

    # using Plots 
    # plot(loading_θs, real.(fp_loading))

    bc1_forward = TractionBoundary(inner=true)
    bc2_forward = TractionBoundary(outer=true)

    loading_profile = BoundaryData(bc1_forward, 
        θs = loading_θs, 
        fields = hcat(fp_loading,fs_loading)
    )

    bd1_for = BoundaryData(ω, bearing, loading_profile)

    basis_order = Int((size(bd1_for.fourier_modes,1) - 1) / 2)
    forward_θs = LinRange(0.0, 2pi, 2*basis_order+2)[1:end-1]

    bd2_for = BoundaryData(bc2_forward, 
        θs = forward_θs, 
        fields = [zeros(Complex{Float64}, 2basis_order + 1) zeros(Complex{Float64}, 2basis_order + 1)]
    )

    modal_method = ModalMethod(tol = 1e-6, only_stable_modes = true)
    forward_sim = BearingSimulation(ω, modal_method, bearing, bd1_for, bd2_for);

    wave = ElasticWave(forward_sim);

    # We should check whether the solution has converged. That is, if the coefficients are getting very small as the fourier mode increase
    l = length(wave.potentials[1].coefficients[1,:])/2 - 1/2
    plot(-l:l,abs.(wave.potentials[1].coefficients[1,:]), xlims = (20,30))

    # Need the frequency mode below to measure the mode 0 of the loading
    frequency_order * bearing.number_of_rollers

# Get the boundary data for the inverse problem from the forward problem
    bc1_inverse = DisplacementBoundary(outer=true)
    bc2_inverse = TractionBoundary(outer=true)

    numberofsensors = 4

    θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]

    # create the data from evaluating the forward problem 
    bd1_inverse = BoundaryData(bc1_inverse, bearing.outer_radius, θs_inv, wave)

    # a little bit of an inverse crime
    bd2_inverse = bd2_for

# Create a fourier basis for the loading, and then create a boundary basis from this

    mZ = frequency_order * bearing.number_of_rollers
    inds = -wave.method.basis_order:wave.method.basis_order
    min_loading_order = minimum(abs.(inds .- mZ))
    max_loading_order = wave.method.basis_order .- mZ

    l = (size(loading_profile.fourier_modes,1) - 1)/2
    loading_profile = fields_to_fouriermodes(loading_profile)
    plot(-l:l,abs.(loading_profile.fourier_modes), xlims = (-10.,10))

    loading_basis_order = 2;
    loading_datas = map(0:loading_basis_order) do n
        fs = zeros(Complex{Float64},2loading_basis_order + 1)
        fs[loading_basis_order + 1 + n] = fs[loading_basis_order + 1 + n] + 0.5
        fs[loading_basis_order + 1 - n] = fs[loading_basis_order + 1 - n] + 0.5

        # The representation of the loading itself
        loading = BoundaryData(bc1_forward, fourier_modes = [fs 0.0 .* fs])

        # the boundary data after considering that the forcing is produced through the contact of the bearings
        BoundaryData(ω, bearing, loading)
    end

    boundarybasis1 = BoundaryBasis(loading_datas)

# solve the inverse problem with the PriorMethod
   method = PriorMethod(tol = modal_method.tol, modal_method = modal_method)

   inverse_sim = BearingSimulation(ω, method, bearing, bd1_inverse, bd2_inverse;
       boundarybasis1 = boundarybasis1,
   );

   inverse_wave = ElasticWave(inverse_sim);

   bd1_inner, bd2_outer = boundary_data(forward_sim, inverse_wave);

   using Plots 
   plot(real.(bd1_inner.fields))
   plot(real.(bd2_outer.fields))

   bd1_inner.fields - bd1_for.fields
end