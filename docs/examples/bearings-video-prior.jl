using ElasticWaves

using Test, Statistics, LinearAlgebra, MultipleScattering
using Plots 

## Settings

    # video and image settings 
    video_width = 420
    image_width = 420

    # the higher the frequency, the worse the result. This is already a high frequency.
    # medium = Elastic(2; ρ = 2.0, cp = 100.0 - 0.001im, cs = 80.0 - 0.001im)
    medium = Elastic(2; ρ = 7000.0, cp = 5000.0, cs = 3500.0)

    # medium = Elastic(2; ρ = 7.0, cp = 4.0 - 0.001im, cs = 3.0 - 0.001im)

    # Ω = 0.25 # the angular speed is normally much smaller than the wavespeeds. But having lower wave speeds makes for nice pictures.

    # Ω = 2pi * 15 / 60 # large wind turbines rotate at about 15 rpm
    Ω = 2pi * 2000 / 60 # large wind turbines rotate at about 15 rpm
    # Ω = 2pi * 3200 / 60 # large wind turbines rotate at about 15 rpm

    Ω * 2.5

    # spread of contact force
    σ = 0.12;

    bearing = RollerBearing(medium = medium, 
        inner_radius = 2.5, outer_radius = 3.5, 
        angular_speed = Ω,  
        rollers_inside = true,
        number_of_rollers = 10,
        roller_radius = 0.5
        # number_of_rollers = 1
    )

    plot(bearing)

    Ω * bearing.outer_radius

    dr = bearing.outer_radius - bearing.inner_radius;

    # Calculate the frequencies to be used for the time video

    # we want enough resolution in time to see the rollers going around. 
    total_time = 2pi / (bearing.number_of_rollers * bearing.angular_speed)

    # to get the number of frames per cycle:
    frames = 12
    dt = total_time / frames

    # from dt we can calculate what's the maximum frequency we need, and then the range of natural frequencies we should use
    maxω = 2pi / dt
    frequency_order = Int(round(maxω / (bearing.number_of_rollers * bearing.angular_speed)))

    frequency_order = 6

    ωms = natural_frequencies(bearing, frequency_order) |> collect

    # need to exclude the zero frequency as is ill posed
    ms = 1:frequency_order
    ωs = ωms
    # ωs = ωms[4:end]

    bc1_forward = TractionBoundary(inner=true)
    bc2_forward = TractionBoundary(outer=true)

    bc1_inverse = DisplacementBoundary(outer=true)
    bc2_inverse = TractionBoundary(outer=true)

## Do one example
    frequency_order = length(ωs) 
    frequency_order = 2
    ω = ωs[frequency_order]

    kp = (ω / medium.cp)
    kp * dr

    loading_basis_order = 15;    
    loading_θs = LinRange(0.0, 2pi, 2loading_basis_order+2)[1:end-1]

    # loading profile
    θo = 3pi/2;
    fp_loading = 0.3 .* exp.(-3.5 .* (sin.(loading_θs) .- sin(θo)).^2) + loading_θs .* 0im; 
 
    # Stribeck
    ε = 0.5;
    fp_loading = (1.0 .- (1.0 + 0.0im .- cos.(loading_θs .- θo)) ./ (2ε)).^(10.0/9.0)
    inds = findall(abs.(imag.(fp_loading)) .> 0.0)
    fp_loading[inds] .= 0.0
    fp_loading = real.(fp_loading)

    fs_loading = 0.0 .* fp_loading;
    

    # add a crack!
    # θo = pi;
    # fp_loading =  fp_loading .- exp.(-20 .* ((loading_θs) .- (θo)).^2) + loading_θs .* 0im; 

    plot(loading_θs, real.(fp_loading))

    θs_loading = (θo - 1.3):0.1:(θo + 1.3)
    scatter!(θs_loading, θs_loading .* 0.0)

    loading_profile = BoundaryData(bc1_forward, θs = loading_θs, 
        fields = hcat(fp_loading,fs_loading)
    )

    loading_modes = -3:3
    loading_modes = -5:5
    loading_modes = -4:4
    loading_profile = fields_to_fouriermodes(loading_profile, loading_modes)
    loading_profile = fouriermodes_to_fields(loading_profile)

    plot(loading_θs, real.(fp_loading))
    plot!(loading_θs, real.(loading_profile.fields[:,1]))

    basis_order = 50;
    θs = LinRange(0.0, 2pi, 2basis_order+2)[1:end-1]

    bd1_for = BoundaryData(ω, bearing, loading_profile; σ = σ)
    bd2_for = BoundaryData(bc2_forward, 
        θs = θs,
        modes = bd1_for.modes,
        coefficients = zeros(Complex{Float64},length(bd1_for.modes),2)
    )
    bd2_for = fouriermodes_to_fields(bd2_for);

    modal_method = ModalMethod(tol = 1e-7, only_stable_modes = true)
    forward_sim = BearingSimulation(ω, bearing, bd1_for, bd2_for; 
        method = modal_method,
        nondimensionalise = true);

    wave = ElasticWave(forward_sim);
    wave.method.mode_errors |> findmax
    
    bd1_inverse = BoundaryData(bc1_inverse, bearing.outer_radius, θs, wave)  
    inverse_fields2 = bd1_inverse.fields[:,1]
    # inverse_fields1 = bd1_inverse.fields[:,1]
    # plot(θs,[real.(inverse_fields2)])
    # plot(θs,[abs.(inverse_fields1),abs.(inverse_fields2)])
    # norm(inverse_fields1 - inverse_fields2) / norm(inverse_fields1)
    # norm(abs.(inverse_fields1) - abs.(inverse_fields2)) / norm(inverse_fields1)


    # Need the frequency mode below to measure the mode 0 of the loading
    frequency_order * bearing.number_of_rollers

    # result_traction = field(wave, bearing, TractionType(); res = 150);
    result_traction = field(wave, bearing, DisplacementType(); res = 150);
    fs = [f[1] for f in result_traction.field];
    result = FrequencySimulationResult(fs, result_traction.x, [ω]);

    field_apply = real
    maxc = 0.5 .* maximum(field_apply.(field(result)))
    minc = min(0.0,field_apply(- maxc))

    gr(size = (550,350))
    plot(result,ω;
        seriestype = :heatmap,
        # legend = :true,
        # phase_time = ts[30],
        # c = :lajolla,
        field_apply = field_apply,
        clims = (minc, maxc),
        # leg = false,
        title = "",
        frame = :none
    )
    mean(norm.(result.field))

    θs_loading = (θo - 1.3):0.1:(θo + 1.3)

    plot!(bearing.outer_radius .* cos.(θs_loading), bearing.outer_radius .* sin.(θs_loading),
        linewidth = 2.0, linestyle = :dash)
    plot!(bearing)
    # savefig("docs/images//bearing-diffraction.png")

## Calculate results for all frequencies

    sims = map(eachindex(ωs)) do i
        println("ω: ", ωs[i])
        println("central basis order: ", (i) * bearing.number_of_rollers)

        bd1_for = BoundaryData(ωs[i], bearing, loading_profile; σ = σ)
        bd2_for = BoundaryData(bc2_forward,
            modes = bd1_for.modes,
            coefficients = zeros(Complex{Float64}, length(bd1_for.modes),2)
        )

        BearingSimulation(ωs[i], bearing, bd1_for, bd2_for; 
            method = modal_method,
            nondimensionalise = true
        )
    end;

    waves = ElasticWave.(sims);

    # check the loading profile
    Z = bearing.number_of_rollers

    loadings = map(eachindex(waves)) do i 
        bd1_inner, bd2_outer = boundary_data(sims[i], waves[i]);
        exp(pi * σ^2 * ms[i]^2) * 2pi / Z .* abs.(bd1_inner.fields[:,1])
    end

    plot(abs.(fp_loading))
    plot!.(loadings)
    plot!()

## Plot results
    results = map(eachindex(ωs)) do i
        result_traction = field(waves[i], bearing, TractionType(); res = 250);
        fs = [f[1] for f in result_traction.field];
        result_pressure = FrequencySimulationResult(fs, result_traction.x, [ωs[i]]);
    end

    results_displace = map(eachindex(ωs)) do i
        result_traction = field(waves[i], bearing, DisplacementType(); res = 250);
        fs = [f[1] for f in result_traction.field];
        result_pressure = FrequencySimulationResult(fs, result_traction.x, [ωs[i]]);
    end

    fields = [[f[1] for f in field(r)] for r in results];
    all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs);

    fields = [[f[1] for f in field(r)] for r in results_displace];
    all_results_displace = FrequencySimulationResult(hcat(fields...), results_displace[1].x, ωs);

    gr(size = (300 * 1.6,300))
    i = 3;
    plot(all_results_displace, ωs[i];
        seriestype=:heatmap
        , field_apply = f -> real(f[1])
    )
    mean(norm.(all_results_displace.field[:,i]))

## Plot a time video
    ts = ω_to_t(ωs)
    # ts = LinRange(0,ts[end],60)
    time_result = frequency_to_time(all_results; t_vec = ts)
    maxc = 0.25 .* maximum(norm.(field(time_result)))
    # maxc = 0.70 .* maximum(norm.(field(time_result)))
    minc = - maxc

    maxc_true_traction = maxc

    t = ts[1]
    
    gr(size = (image_width, image_width))
    plot(time_result, t,
        seriestype=:heatmap,
        # field_apply = f -> real(f[2]),
        clim = (minc, maxc),
        leg = false,
    )
    # scatter!([bearing.outer_radius * cos(θo)], [bearing.outer_radius * sin(θo)])
    plot!(bearing, t)
    plot!(frame = :none, title="", xguide ="",yguide ="")

    # savefig("docs/images/bearings-true.pdf")


    gr(size = (video_width, video_width))
    anim = @animate for t in ts[1:end]
        plot(time_result, t,
            seriestype=:heatmap,
            # field_apply = f -> real(f[2]),
            clim = (minc, maxc),
            leg = false,
        )
        # scatter!([bearing.outer_radius * cos(θo)], [bearing.outer_radius * sin(θo)])
        plot!(bearing, t)
        plot!(frame = :none, title="", xguide ="",yguide ="")
    end
    # gif(anim,"bearings-true-dynamic.gif", fps = 4)
    # gif(anim,"docs/images/bearings-true.gif", fps = 4)
    
    time_result = frequency_to_time(all_results_displace; t_vec = ts)
    maxc = 0.70 .* maximum(norm.(field(time_result)))
    # maxc = 0.70 .* maximum(norm.(field(time_result)))
    minc = - maxc

    maxc_true_displace = maxc

    t = ts[1]
    
    gr(size = (image_width, image_width))
    plot(time_result, t,
        seriestype=:heatmap,
        # field_apply = f -> real(f[2]),
        clim = (minc, maxc),
        leg = false,
    )
    # scatter!([bearing.outer_radius * cos(θo)], [bearing.outer_radius * sin(θo)])
    plot!(bearing, t)
    plot!(frame = :none, title="", xguide ="",yguide ="")

    # savefig("docs/images/bearings-true-displace.pdf")

    gr(size = (300 * 1.6,300))

    anim = @animate for t in ts[1:end]
        plot(time_result, t,
            seriestype=:heatmap,
            # field_apply = f -> real(f[2]),
            clim = (minc, maxc),
            leg = false,
        )
        # scatter!([bearing.outer_radius * cos(θo)], [bearing.outer_radius * sin(θo)])
        plot!(bearing, t)
        plot!(frame = :none, title="", xguide ="",yguide ="")
    end
    # gif(anim,"bearings-true-dynamic-displace.gif", fps = 4)
    # gif(anim,"docs/images/bearings-true-displace.gif", fps = 4)

## Use loading prior 

    i = 1;
    i = length(ωs);
    ω = ωs[i]

    # For steel, this problem is quite sensitive to errors and changes in the boundary data. Naturally, a statistical method is needed here to regularise or complete avoid unstable problems. 
    
    numberofsensors = 4
    # numberofsensors = 8
    # numberofsensors = 12
    # numberofsensors = 16

    θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]
    # θs_inv = [3.0,4.0,5.0]

    # create the data from evaluating the forward problem
    bd1_inverse = BoundaryData(bc1_inverse, bearing.outer_radius, θs_inv, waves[i])

    # bd2_inverse = BoundaryData(bc2_inverse, 
    #     θs = θs,
    #     fields = zeros(Complex{Float64},length(θs),2)
    # )

    # Create a fourier basis for the loading, and then create a boundary basis from this

    # as min_mode = 0 and max_mode is very large, we can consider any loading profiles with modes = -max_modes:max_modes 

    loading_basis_order = 2;
    basis = map(-loading_basis_order:loading_basis_order) do n
        fp = [1.0 + 0.0im]

        # The representation of the loading itself
        loading_data = BoundaryData(bc1_forward, 
            modes = [n],
            coefficients = [fp 0.0.*fp]
        )

        # How loading is translated into boundary data
        BoundaryData(ω, bearing, loading_data; σ = σ)
    end;
    boundarybasis1 = BoundaryBasis(basis);

    bd2_inverse = BoundaryData(bc2_inverse, 
        modes = boundarybasis1.basis[1].modes,
        coefficients = zeros(Complex{Float64},length(boundarybasis1.basis[1].modes),2)
    )

    # solve the inverse problem with the PriorMethod
    method = PriorMethod(tol = modal_method.tol, modal_method = modal_method)
    
    inverse_sim = BearingSimulation(ω, method, bearing, bd1_inverse, bd2_inverse;
        boundarybasis1 = boundarybasis1
    );

    inverse_wave = ElasticWave(inverse_sim);
    inverse_wave.method.condition_number
    inverse_wave.method.boundary_error

    inv_modes = inverse_wave.potentials[1].modes

    inds = [findfirst(m .== waves[i].potentials[1].modes) for m in inv_modes]
    
    waves[i].potentials[1].modes[inds] == inv_modes

    bd1_inner, bd2_outer = boundary_data(sims[i], inverse_wave);
    inv_loading = exp(pi*σ^2*ms[i]^2) * 2pi / Z .* abs.(bd1_inner.fields[:,1])

    gr(size = (image_width * 1.6, image_width))
    plot(loading_θs, abs.(loadings[1]), label = "Original")
    plot!(loading_θs, inv_loading, label = "inv. with prior", linestyle = :dash) 
    plot!(ylab = "loading profile", xlab = "θ")
    
    # wave = waves[i];
    # sim = inverse_sim

    errors = [
        norm(waves[i].potentials[1].coefficients[:,inds[m]] - inverse_wave.potentials[1].coefficients[:,m]) / norm(waves[i].potentials[1].coefficients[:,inds[m]])
    for m = eachindex(inds)];

    maximum(errors)
    
    ## do all frequencies
    inverse_waves = map(eachindex(ωs)) do i

        println("ω: ", ωs[i])
        # println("central basis order: ", (i+1) * bearing.number_of_rollers)

         # create the data from evaluating the forward problem
        bd1_inverse = BoundaryData(bc1_inverse, bearing.outer_radius, θs_inv, waves[i])

        basis = map(-loading_basis_order:loading_basis_order) do n
            fp = [1.0 + 0.0im]
    
            # The representation of the loading itself
            loading_data = BoundaryData(bc1_forward, 
                modes = [n],
                coefficients = [fp 0.0.*fp]
            )
    
            # How loading is translated into boundary data
            BoundaryData(ωs[i], bearing, loading_data; σ = σ)
        end;

        boundarybasis1 = BoundaryBasis(basis);
    
        bd2_inverse = BoundaryData(bc2_inverse, 
            modes = boundarybasis1.basis[1].modes,
            coefficients = zeros(Complex{Float64},length(boundarybasis1.basis[1].modes),2)
        )
    
        # solve the inverse problem with the PriorMethod
        method = PriorMethod(tol = modal_method.tol, modal_method = modal_method)
    
        inverse_sim = BearingSimulation(ωs[i], method, bearing, bd1_inverse, bd2_inverse;
            boundarybasis1 = boundarybasis1
        );

        inverse_wave = ElasticWave(inverse_sim)

        println("Condition num: ", inverse_wave.method.condition_number)
        println("Boundary error: ", inverse_wave.method.boundary_error)

        inverse_wave
    end;

    inv_loadings = map(eachindex(inverse_waves)) do i 
        bd1_inner, bd2_outer = boundary_data(sims[i], inverse_waves[i]);
        exp(pi*σ^2*ms[i]^2) * 2pi / Z .* abs.(bd1_inner.fields[:,1])
    end;

    gr(size = (image_width * 1.6, image_width))
    plot(loading_θs, abs.(loadings[1]), label = "Original")
    [
        plot!(loading_θs, l, label = "inv. with prior", linestyle = :dash) 
    for l in inv_loadings]
    plot!(ylab = "loading profile", xlab = "θ")

    gr(size = (image_width, image_width / 1.6))
    plot(loading_θs, abs.(loadings[1]), label = "Original")
    plot!(loading_θs, mean(inv_loadings), label = "inv. with prior", linestyle = :dash) 

    # savefig("load-profile-inv-with-prior.pdf")

## Plot results
    results = map(eachindex(ωs)) do i
        result_traction = field(inverse_waves[i], bearing, TractionType(); res = 250);
        fs = [f[1] for f in result_traction.field];
        result_pressure = FrequencySimulationResult(fs, result_traction.x, [ωs[i]]);
    end

    results_displace = map(eachindex(ωs)) do i
        result_traction = field(inverse_waves[i], bearing, DisplacementType(); res = 250);
        fs = [f[1] for f in result_traction.field];
        result_pressure = FrequencySimulationResult(fs, result_traction.x, [ωs[i]]);
    end

    fields = [[f[1] for f in field(r)] for r in results];
    all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs);

    fields = [[f[1] for f in field(r)] for r in results_displace];
    all_results_displace = FrequencySimulationResult(hcat(fields...), results_displace[1].x, ωs);

    gr(size = (350,300))
    # i = 10;
    plot(all_results, ωs[end-3];
        seriestype=:heatmap
        , field_apply = f -> real(f[1])
    )

## Plot a time video
    ts = ω_to_t(ωs)

    time_result = frequency_to_time(all_results; t_vec = ts)
    maxc = 0.25 .* maximum(norm.(field(time_result)))

    maxc_true_traction
    maxc = maxc_true_traction
    minc = - maxc

    t = ts[1]
    gr(size = (image_width, image_width))
    plot(time_result, t,
        seriestype=:heatmap,
        # field_apply = f -> real(f[2]),
        clim = (minc, maxc),
        leg = false,
    )
    scatter!([bearing.outer_radius .* cos.(θs_inv)], [bearing.outer_radius .* sin.(θs_inv)])
    plot!(bearing, t)
    plot!(frame = :none, title="", xguide ="",yguide ="")
    # savefig("docs/images/bearings-with-prior-traction.pdf")
    
    gr(size = (video_width, video_width))


    anim = @animate for t in ts[1:end]
        plot(time_result, t,
            seriestype=:heatmap,
            # field_apply = f -> real(f[2]),
            clim = (minc, maxc),
            leg = false,
        )
        scatter!([bearing.outer_radius .* cos.(θs_inv)], [bearing.outer_radius .* sin.(θs_inv)])
        plot!(bearing, t)
        plot!(frame = :none, title="", xguide ="",yguide ="")
    end

    # gif(anim,"docs/images/bearings-with-prior-traction.gif", fps = 4)

    time_result = frequency_to_time(all_results_displace; t_vec = ts)
    maxc = maxc_true_displace
    minc = - maxc

    gr(size = (image_width, image_width))
    plot(time_result, t,
        seriestype=:heatmap,
        # field_apply = f -> real(f[1]),
        clim = (minc, maxc),
        leg = false,
    )
    scatter!([bearing.outer_radius .* cos.(θs_inv)], [bearing.outer_radius .* sin.(θs_inv)])
    plot!(bearing, t)
    plot!(frame = :none, title="", xguide ="",yguide ="")
    # savefig("docs/images/bearings-with-prior-displacement.pdf")

    t = ts[4]
    anim = @animate for t in ts[1:end]
        plot(time_result, t,
            seriestype=:heatmap,
            # field_apply = f -> real(f[2]),
            clim = (minc, maxc),
            leg = false,
        )
        scatter!([bearing.outer_radius .* cos.(θs_inv)], [bearing.outer_radius .* sin.(θs_inv)])
        plot!(bearing, t)
        plot!(frame = :none, title="", xguide ="",yguide ="")
    end
    # gif(anim,"docs/images/bearings-with-prior-displace.gif", fps = 4)