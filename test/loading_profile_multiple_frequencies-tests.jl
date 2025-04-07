using ElasticWaves, Test, LinearAlgebra

# the higher the frequency, the worse the result. This is already a high frequency.
medium = Elastic(2; ρ = 2.0, cp = 10.0 - 0.0im, cs = 8.0 - 0.0im)

Ω = 2pi * 15 / 60 # large wind turbines rotate at about 15 rpm
Z = 8 

bearing = RollerBearing(medium = medium, 
    inner_radius = 1.5, outer_radius = 2.0, 
    angular_speed = Ω,  
    rollers_inside = true,
    number_of_rollers = Z
)

frequency_order = 4

ωms = natural_frequencies(bearing, frequency_order) |> collect

ωms = [ωms[4]]
ω = ωms[1]

dr = bearing.outer_radius - bearing.inner_radius
kp = (ω / medium.cp)
kp * dr

# create the true loading profile, then solve the forward problem to create data for the inverse problem
    bc1_forward = TractionBoundary(inner=true)
    bc2_forward = TractionBoundary(outer=true)

    loading_basis_order_forward = 2; # Used for the forward problem
    loading_basis_order = 2; # Used for the inverse problem

    # loading_profile2 = fields_to_fouriermodes(loading_profile, -loading_basis_order_forward:loading_basis_order_forward)
    # loading_profile = fouriermodes_to_fields(loading_profile2)

    loading_profile = BoundaryData(bc1_forward,
        modes = [-2,-1,0,1,2],
        coefficients = [
            0.3 0.0im;
            -0.2 0.0im;
            0.0 0.0im;
            -0.2 0.0im;
            0.3 0.0im
        ]
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
    # println("max mode error:",maximum(waves[1].method.mode_errors .|> abs))

    # Need the frequency mode below to measure the mode 0 of the loading
    frequency_order * bearing.number_of_rollers

# Get the boundary data for the inverse problem from the forward problem
    bc1_inverse = DisplacementBoundary(outer=true)
    bc2_inverse = TractionBoundary(outer=true)

    numberofsensors = 4

    θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]

    loading_resolution = 120;
    loading_θs = LinRange(0.0, 2pi, 2*loading_resolution+2)[1:end-1]

    # solve the inverse problem with the PriorMethod
    method = PriorMethod(tol = modal_method.tol, modal_method = modal_method)

    inverse_sims = map(eachindex(ωms)) do i

        ω = ωms[i]

        # Create a fourier basis for the loading, and then create a boundary basis from this
        loading_basis_order = 2;  
        basis = map(-loading_basis_order:loading_basis_order) do n
            fp = [1.0 + 0.0im]

            # The representation of the loading itself
            loading_data = BoundaryData(bc1_forward, 
                modes = [n],
                coefficients = [fp 0.0.*fp]
            )

            # How loading is translated into boundary data
            BoundaryData(ω, bearing, loading_data)
        end;
        boundarybasis1 = BoundaryBasis(basis);

        # create the data from evaluating the forward problem 
        bd1_inverse = BoundaryData(bc1_inverse, bearing.outer_radius, θs_inv, waves[i])

        bd2_inverse = BoundaryData(bc2_inverse, 
            modes = boundarybasis1.basis[1].modes,
            coefficients = 0.0 .* boundarybasis1.basis[1].coefficients
        )
        
        inverse_sim = BearingSimulation(ω, method, bearing, bd1_inverse, bd2_inverse;
            boundarybasis1 = boundarybasis1,
        )
    end;

    inverse_wave = ElasticWave(inverse_sims[1]);
    inverse_wave.method.condition_number
    inverse_wave.method.boundary_error

    bd1_inner, bd2_outer = boundary_data(forward_sims[1], inverse_wave);
    norm(bd1_inner.fields - forward_sims[1].boundarydata1.fields) / norm(forward_sims[1].boundarydata1.fields)

    inverse_sim_copies = deepcopy.(inverse_sims);
    nondimensionalise!.(inverse_sim_copies);

    import ElasticWaves: prior_and_bias_inverse
    
    data = prior_and_bias_inverse.(inverse_sim_copies);

    E_invs = [d[1] for d in data];
    EBs = [d[1] * d[2] for d in data];
    Eds = [d[1] * d[3] for d in data];
    y_invs = [d[4] for d in data];

    BB = vcat(EBs...);
    DD = vcat(Eds...);
    YY = vcat(y_invs...);

    cond(BB)

    # Did non-dimenalisation make x different for each frequency??
    x = BB \ (YY - DD);

    # This should re-dimenalisation x as it has units of traction
    λ2μ = bearing.medium.ρ * bearing.medium.cp^2

    -loading_basis_order:loading_basis_order
    x = x .* λ2μ


    loading_profile2.modes
    loading_profile2.coefficients

    # δ = method.regularisation_parameter
    # bigA = [BB; sqrt(δ) * I];
    # x = bigA \ [YY - DD; zeros(size(BB)[2])]


    # using Plots
    # scatter(wave.potentials[1].modes, abs.(wave.potentials[1].coefficients[1,:]), markersize = 4.0)
    # scatter!(inverse_wave.potentials[1].modes, abs.(inverse_wave.potentials[1].coefficients[1,:]), markersize = 3.0)
    # plot!(xlims = (-1 -loading_basis_order - mZ, loading_basis_order - mZ + 1))

    # predicted_forcing_coefficients = hcat(
    #     field_modes(inverse_wave, bearing.inner_radius, bd1_for.boundarytype.fieldtype),
    #     field_modes(inverse_wave, bearing.outer_radius, bd1_for.boundarytype.fieldtype)
    # )

    # scatter(inverse_wave.method.modal_method.modes, abs.(predicted_forcing_coefficients[:,1]))
    # scatter!(inverse_wave.method.modal_method.modes, abs.(predicted_forcing_coefficients[:,1]))

#    plot(bd1_inner.θs, 2pi / Z .* abs.(bd1_inner.fields[:,1]), label = "predicted loading")
#    plot!(loading_θs, abs.(fp_loading), linestyle = :dash, label = "true loading")
   
#    plot(real.(bd2_outer.fields))