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

frequency_order = 3

ωms = natural_frequencies(bearing, frequency_order) |> collect

ωms = [ωms[frequency_order]]
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

    # loading_resolution = 120;
    # loading_θs = LinRange(0.0, 2pi, 2*loading_resolution+2)[1:end-1]

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
        bd1_inverse = BoundaryData(
            BoundaryData(bc1_inverse; θs = θs_inv), 
            bearing.outer_radius, waves[i]
        )

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
    norm(bd1_inner.modes - forward_sims[1].boundarydata1.modes) 
    norm(bd2_outer.modes - forward_sims[1].boundarydata2.modes)

    norm(bd1_inner.coefficients - forward_sims[1].boundarydata1.coefficients) / norm(forward_sims[1].boundarydata1.coefficients)
    
    norm(bd2_outer.coefficients - forward_sims[1].boundarydata2.coefficients) 
    
    σ = bearing.roller_contact_angular_spread;
    scale = Z/(2pi) * exp(- (π * σ^2) *  frequency_order^2);
    
    loading_ordered = select_modes(loading_profile, bd1_inner.modes .- frequency_order * Z);

    norm(loading_ordered.coefficients - bd1_inner.coefficients ./ scale) / norm(loading_ordered.coefficients)
    bd1_inner.modes .- frequency_order * Z

    inverse_sim_copies = deepcopy.(inverse_sims);
    nondimensionalise!.(inverse_sim_copies);

    import ElasticWaves: prior_and_bias_inverse
    
    data = prior_and_bias_inverse(inverse_sims[1]);

    E_inv = data[1];
    B = data[2];
    d = data[3];
    y_inv = data[4];

    data = prior_and_bias_inverse(inverse_sim_copies[1]);

    E_inv_non = data[1];
    B_non = data[2];
    d_non = data[3];
    y_inv_non = data[4];

    kP = ω / bearing.medium.cp;
    kS = ω / bearing.medium.cs;

    ρλ2μ = bearing.medium.ρ * bearing.medium.cp^2

    norm(B * kP^2 * ρλ2μ - B_non) / norm(B_non)
    norm(y_inv .* kp - y_inv_non) / norm(y_inv_non)
    norm(d - d_non .* ρλ2μ)


    data = prior_and_bias_inverse.(inverse_sim_copies);

    E_invs = [d[1] for d in data];
    Bs = [d[2] for d in data];
    ds = [d[3] for d in data];

    EBs = [d[1] * d[2] for d in data];
    Eds = [d[1] * d[3] for d in data];
    y_invs = [d[4] for d in data];

    BB = vcat(EBs...);
    DD = vcat(Eds...);
    YY = vcat(y_invs...);

    cond(BB)

    # Did non-dimenalisation make x different for each frequency??
    x = BB \ (YY - DD);
    modes = inverse_sim_copies[1].method.modal_method.modes

    coefficients = Bs[1] * x + ds[1];
    coefficients = reshape(coefficients,(4,:)) |> transpose |> collect;

    kP = ω / bearing.medium.cp;
    kS = ω / bearing.medium.cs;

    coefficients = coefficients ./ kP^2

    pressure_coefficients = coefficients[:,1:2] |> transpose
    shear_coefficients = coefficients[:,3:4] |> transpose

    φ = HelmholtzPotential(bearing.medium.cp, kP, pressure_coefficients, modes)
    ψ = HelmholtzPotential(bearing.medium.cs, kS, shear_coefficients, modes)

    inverse_wave = ElasticWave(ω, bearing.medium, [φ, ψ], method);

    bd1_inner, bd2_outer = boundary_data(forward_sims[1], inverse_wave);
    
    norm(bd1_inner.coefficients - forward_sims[1].boundarydata1.coefficients) / norm(forward_sims[1].boundarydata1.coefficients)

    # This should re-dimenalisation x as it has units of traction
    λ2μ = bearing.medium.ρ * bearing.medium.cp^2

    loading_profile2 = select_modes(loading_profile, -loading_basis_order:loading_basis_order)
    loading_profile2.modes
    norm(x .* λ2μ - loading_profile2.coefficients[:,1])

