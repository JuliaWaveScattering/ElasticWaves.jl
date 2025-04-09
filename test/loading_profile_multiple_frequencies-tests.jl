# using ElasticWaves, Test, LinearAlgebra, Statistics

@testset "Loading profile multiple frequencies" begin

# the higher the frequency, the worse the result. This is already a high frequency.
medium = Elastic(2; ρ = 2.0, cp = 10.0 - 0.0im, cs = 8.0 - 0.0im)

Ω = 2pi * 15 / 60 # large wind turbines rotate at about 15 rpm
Ω = 2pi * 60 / 60 # large wind turbines rotate at about 15 rpm
Z = 8 

ratio_shear_to_normal = 0.0

bearing = RollerBearing(medium = medium, 
    inner_radius = 1.5, outer_radius = 2.0, 
    angular_speed = Ω,  
    rollers_inside = true,
    number_of_rollers = Z
)

frequency_order = 12

ωms = natural_frequencies(bearing, frequency_order) |> collect

# ωms = [ωms[frequency_order]]
ω = ωms[1]

dr = bearing.outer_radius - bearing.inner_radius
kp = (ω / medium.cp)
kp * dr

# create the true loading profile, then solve the forward problem to create data for the inverse problem
    bc1_forward = TractionBoundary(inner=true)
    bc2_forward = TractionBoundary(outer=true)

    loading_basis_order_forward = 2; # Used for the forward problem
    loading_basis_order = 2; # Used for the inverse problem

    loading_coefficients = [0.3, -0.2, 0.0, -0.2, 0.3]
    loading_profile = BoundaryData(bc1_forward,
        modes = [-2,-1,0,1,2],
        coefficients = hcat(loading_coefficients, ratio_shear_to_normal .* loading_coefficients)
    )

    bd1_for = BoundaryData(ω, bearing, loading_profile)

    modal_method = ModalMethod(tol = 1e-9, only_stable_modes = true)
    forward_sims = map(ωms) do ω
        bd1_for = BoundaryData(ω, bearing, loading_profile);
        bd2_for = BoundaryData(bc2_forward, 
            modes = bd1_for.modes, 
            coefficients = 0.0 .* bd1_for.coefficients
        )

        forward_sim = BearingSimulation(ω, modal_method, bearing, bd1_for, bd2_for);
    end;

    waves = ElasticWave.(forward_sims);
    println("max mode error:",maximum(waves[1].method.mode_errors .|> abs))
    println("max mode error:",maximum(waves[end].method.mode_errors .|> abs))

    # Need the frequency mode below to measure the mode 0 of the loading
    frequency_order * bearing.number_of_rollers

# Get the boundary data for the inverse problem from the forward problem
    bc1_inverse = DisplacementBoundary(outer=true)
    bc2_inverse = TractionBoundary(outer=true)

    numberofsensors = 1

    θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]

    method = ConstantRollerSpeedMethod(
        tol = modal_method.tol, modal_method = modal_method,
        loading_modes = -2:2,
        ratio_shear_to_normal = ratio_shear_to_normal,
    )

    function inverse_simulations(error = 0.0)
        inverse_sims = map(eachindex(ωms)) do i

            # Create a fourier basis for the loading, and then create a boundary basis from this
            boundarybasis1 = BoundaryBasis(ωms[i], bearing, method);

            # create boundary data by evaluating the forward problem 
            bd1_inverse = BoundaryData(
                BoundaryData(bc1_inverse; θs = θs_inv), 
                bearing.outer_radius, waves[i]
            )

            #add error to the measurements
            scale = error * mean(abs.(bd1_inverse.fields)) / 2.0;
            field_error = scale .* (rand(size(bd1_inverse.fields)...) .- 0.5) + scale .* (rand(size(bd1_inverse.fields)...) .- 0.5) .* im

            bd1_inverse = BoundaryData(bc1_inverse,
                θs = bd1_inverse.θs,
                fields = bd1_inverse.fields + field_error
            )

            # known zero traction bvoundary data
            bd2_inverse = BoundaryData(bc2_inverse,
                modes = boundarybasis1.basis[1].modes,
                coefficients = 0.0 .* boundarybasis1.basis[1].coefficients
            )
            
            inverse_sim = BearingSimulation(ωms[i], method, bearing, bd1_inverse, bd2_inverse;
                boundarybasis1 = boundarybasis1,
            )
        end
    end

    # function ElasticWaves(sims::Vector{B}, method::AbstractPriorMethod) where B <: BearingSimulation

    inverse_waves = ElasticWaveVector(inverse_simulations());

    loading_profile2 = select_modes(loading_profile, inverse_waves[1].method.loading_modes)
    
    error = norm(inverse_waves[1].method.loading_coefficients - loading_profile2.coefficients[:,1]) / norm(inverse_waves[1].method.loading_coefficients)

    @test error < 1e-12

    inverse_waves = ElasticWaveVector(inverse_simulations(0.01));

    loading_profile2 = select_modes(loading_profile, inverse_waves[1].method.loading_modes)
    
    error = norm(inverse_waves[1].method.loading_coefficients - loading_profile2.coefficients[:,1]) / norm(inverse_waves[1].method.loading_coefficients)

    @test error < 0.065

    # loading_profile_predict = BoundaryData(bc1_forward,
    #     modes = inverse_waves[1].method.loading_modes,
    #     coefficients = hcat(
    #         inverse_waves[1].method.loading_coefficients, 
    #         ratio_shear_to_normal .* inverse_waves[1].method.loading_coefficients
    #     )
    # )

    # θs = 0.0:0.1:(2pi)
    # loading_profile = fouriermodes_to_fields(loading_profile,θs)
    # loading_profile_predict = fouriermodes_to_fields(loading_profile_predict,θs)

    # using Plots
    # gr(size = (320,240), linewidth = 2.0)
    # plot(bearing)

    # plot(loading_profile.θs, real.(loading_profile.fields[:,1]))
    # plot!(loading_profile_predict.θs, real.(loading_profile_predict.fields[:,1]), linestyle =:dash)
    
 
end