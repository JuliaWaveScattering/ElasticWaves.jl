using ElasticWaves, Test, LinearAlgebra, Statistics

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

    inverse_sims = map(eachindex(ωms)) do i

        # Create a fourier basis for the loading, and then create a boundary basis from this
        boundarybasis1 = BoundaryBasis(ωms[i], bearing, method);

        # create boundary data by evaluating the forward problem 
        bd1_inverse = BoundaryData(
            BoundaryData(bc1_inverse; θs = θs_inv), 
            bearing.outer_radius, waves[i]
        )

        # known zero traction bvoundary data
        bd2_inverse = BoundaryData(bc2_inverse,
            modes = boundarybasis1.basis[1].modes,
            coefficients = 0.0 .* boundarybasis1.basis[1].coefficients
        )
        
        inverse_sim = BearingSimulation(ωms[i], method, bearing, bd1_inverse, bd2_inverse;
            boundarybasis1 = boundarybasis1,
        )
    end;

    # function ElasticWaves(sims::Vector{B}, method::AbstractPriorMethod) where B <: BearingSimulation

    inverse_waves = ElasticWaveVector(inverse_sims, method);


    # This should re-dimenalisation x as it has units of traction
    λ2μ = bearing.medium.ρ * bearing.medium.cp^2

    loading_profile2 = select_modes(loading_profile, -method.loading_basis_order:method.loading_basis_order)
    loading_profile2.modes
    norm(x .* λ2μ - loading_profile2.coefficients[:,1]) / norm(x .* λ2μ)



    inverse_sim_copies = deepcopy.(inverse_sims);
    nondimensionalise!.(inverse_sim_copies);

    import ElasticWaves: prior_and_bias_inverse
    
    data = prior_and_bias_inverse.(inverse_sim_copies);

    Bs = [d[2] for d in data];
    ds = [d[3] for d in data];

    EBs = [d[1] * d[2] for d in data];
    Eds = [d[1] * d[3] for d in data];
    y_invs = [d[4] for d in data];

    # add some error to the boundary data
    map(y_invs) do y 
        y[1:2] = y[1:2] + mean(abs.(y[1:2])) .* 0.001 .* (rand(2) .- 0.5) + mean(abs.(y[1:2])) .* 0.001 .* (rand(2) .- 0.5) .* im
    end;

    BB = vcat(EBs...);
    DD = vcat(Eds...);
    YY = vcat(y_invs...);

    cond(BB)

    SM = diagm([4.0 / sum(abs.(BB[:,j])) for j in 1:size(BB,2)])
    BB = BB * SM

    cond(BB)

    x = BB \ (YY - DD);
    x = SM * x

    if size(BB,1) < size(BB,2)
        @error "For the constant speed roller method, we expected the system to recover the Fourier coefficients of the loading profile to be overdetermined, but it is not. Consider using more frequencies."
    end     

    # This should re-dimenalisation x as it has units of traction
    λ2μ = bearing.medium.ρ * bearing.medium.cp^2

    loading_profile2 = select_modes(loading_profile, -method.loading_basis_order:method.loading_basis_order)
    loading_profile2.modes
    norm(x .* λ2μ - loading_profile2.coefficients[:,1]) / norm(x .* λ2μ)


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



