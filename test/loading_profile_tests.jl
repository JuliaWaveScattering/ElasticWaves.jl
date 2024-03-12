# This is very commented example of how to use the PriorMethod

@testset "Loading profile" begin

# the higher the frequency, the worse the result. This is already a high frequency.
medium = Elastic(2; ρ = 2.0, cp = 1.0 - 0.0im, cs = 0.8 - 0.0im)

Ω = 0.02 # the angular speed is normally much smaller than the wavespeeds. But having lower wave speeds makes for nice pictures.

bearing = RollerBearing(medium = medium, 
    inner_radius = 1.5, outer_radius = 2.0, 
    angular_speed = Ω,  
    rollers_inside = true
)

frequency_order = 4

ωms = natural_frequencies(bearing, frequency_order) |> collect

ω = ωms[end]

dr = bearing.outer_radius - bearing.inner_radius
kp = (ω / medium.cp)
kp * dr

# create the true loading profile, then solve the forward problem to create dat for tthe inverse problem

    loading_resolution = 30;
    loading_θs = LinRange(0.0, 2pi, 2*loading_resolution+2)[1:end-1]

    θo = 3pi/2;
    fp_loading = 0.2 .- exp.(-0.4 .* (sin.(loading_θs) .- sin(θo)).^2) + loading_θs .* 0im; 
    # fp_loading = -0.3 .+ exp.(-1.0 .* (loading_θs .- θo).^2) + loading_θs .* 0im; 
    fs_loading = 0.1 .* fp_loading;

    # using Plots 
    # plot(loading_θs, real.(fp_loading))

    loading_basis_order = 1
    bc1_forward = TractionBoundary(inner=true)
    bc2_forward = TractionBoundary(outer=true)

    loading_profile = BoundaryData(bc1_forward, 
        θs = loading_θs, 
        fields = hcat(fp_loading,fs_loading)
    )

    bd1_for = BoundaryData(ω, bearing, loading_profile)
    bd2_for = BoundaryData(bc2_forward, 
        θs = loading_θs, 
        fourier_modes = 0.0 .* bd1_for.fourier_modes
    )

    modal_method = ModalMethod(tol = 1e-9, only_stable_modes = true)
    forward_sim = BearingSimulation(ω, bearing, bd1_for, bd2_for; 
        method = modal_method,
        nondimensionalise = true);

    wave = ElasticWave(forward_sim);

# Get the boundary data for the inverse problem from the forward problem
    bc1_inverse = DisplacementBoundary(outer=true)
    bc2_inverse = TractionBoundary(outer=true)

    numberofsensors = 3
        
    θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]

    # create the data from evaluating the forward problem 
    bd1_inverse = BoundaryData(bc1_inverse, bearing.outer_radius, θs_inv, wave)

    # a little bit of an inverse crime
    bd2_inverse = bd2_for

# Create a fourier basis for the loading, and then create a boundary basis from this

    
    loading_profile = BoundaryData(bc1_forward, fields = hcat(fp_loading,fs_loading))

    basis_order = 15;
    θs = LinRange(0.0, 2pi, 2basis_order+2)[1:end-1]

    bd1_for = BoundaryData(ω, bearing, loading_profile)
end