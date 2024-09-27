using ElasticWaves

using Test, Statistics, LinearAlgebra, MultipleScattering
using Plots 

# the higher the frequency, the worse the result. This is already a high frequency.
medium = Elastic(2; ρ = 1.0, cp = 1.0 - 0.01im, cs = 0.8 - 0.01im)

# Using realistic value for steel lead to needing to solve the low frequency limit, which currently does not work
# medium = Elastic(2; ρ = 7000.0, cp = 5000.0 - 0.0im, cs = 3500.0 - 0.0im)

Ω = 2pi * 10 / 60 
Ω = 2pi * 0.2 / 60
Z = 12; 

# spread of contact force
σ = 0.2;

bearing = RollerBearing(medium = medium, 
    inner_radius = 10.0, outer_radius = 15.0, 
    angular_speed = Ω,  
    rollers_inside = true,
    number_of_rollers = Z,
    roller_radius = 1.8
)

Ω * bearing.outer_radius

plot(bearing, 0.0)

# Calculate the frequencies to be used for the time video

    # we want enough resolution in time to see the rollers going around. 
    total_time = 2pi / (bearing.angular_speed)

    # to get the number of frames per cycle:
    frames = 30
    dt = total_time / frames

    # from dt we can calculate what's the maximum frequency we need, and then the range of natural frequencies we should use
    maxω = 2pi / dt
    frequency_order = Int(round(maxω / (bearing.number_of_rollers * bearing.angular_speed)))

    # frequency_order = 100

    ms = 1:frequency_order |> collect
    ωms = natural_frequencies(bearing, frequency_order) |> collect
    Z * Ω 

    ωs = ωms

    dr = bearing.outer_radius - bearing.inner_radius
    kps = (ωms ./ medium.cp)
    krs = kps .* bearing.outer_radius

    # inds = findall(abs.(krs) .<= 20.0)
    # ωms = ωms[inds]

## Do one example
    ω = ωms[end-2]
    ω = ωms[1]
    m = ms[1]

    loading_basis_order = 15;    
    loading_θs = LinRange(0.0, 2pi, 2loading_basis_order+2)[1:end-1]

    # loading profile
    θo = 3pi/2;
    fp_loading = 0.3 .* exp.(-2. .* (sin.(loading_θs) .- sin(θo)).^2) + loading_θs .* 0im; 
    # fp_loading = -0.3 .+ exp.(-1.0 .* (loading_θs .- θo).^2) + loading_θs .* 0im; 
    fs_loading = 0.1 .* fp_loading;

    # add a crack!
    # θo = pi;
    # fp_loading =  fp_loading .- exp.(-20 .* ((loading_θs) .- (θo)).^2) + loading_θs .* 0im; 

    plot(loading_θs, real.(fp_loading))

    bc1_forward = TractionBoundary(inner=true)
    bc2_forward = TractionBoundary(outer=true)

    loading_profile = BoundaryData(bc1_forward, θs = loading_θs, fields = hcat(fp_loading,fs_loading));


    bd1_for = BoundaryData(ω, bearing, loading_profile; σ = σ);   

    basis_order = maximum(abs.(bd1_for.modes)) + 5;
    θs = LinRange(0.0, 2pi, 2basis_order+2)[1:end-1]

    bd2_for = BoundaryData(bc2_forward, θs=θs, fields = hcat(0.0.*θs .+ 0.0im, 0.0.*θs .+ 0.0im))

    modal_method = ModalMethod(tol = 4e-4, only_stable_modes = true)
    forward_sim = BearingSimulation(ω, bearing, bd1_for, bd2_for; 
        method = modal_method);

    wave = ElasticWave(forward_sim);
    bd1_for.modes |> mean

    # test if loading profile was recovered
    # bd1_inner, bd2_outer = boundary_data(forward_sim, wave);
    bd1_inner = BoundaryData(bd1_for.boundarytype, bearing.inner_radius, LinRange(0,2pi,200)[1:end-1], wave);
    plot(bd1_inner.θs, exp(pi*σ^2*m^2) * 2pi / bearing.number_of_rollers .* abs.(bd1_inner.fields[:,1]), label = "predicted loading")
    plot!(loading_θs, abs.(fp_loading), linestyle = :dash, label = "true loading")

## plot the whole field    
    result_traction = field(wave, bearing, TractionType(); res = 150);
    fs = [f[1] for f in result_traction.field];
    result_pressure = FrequencySimulationResult(fs, result_traction.x, [ω]);

    field_apply = abs
    
    plot(result_pressure,ω;
        seriestype = :heatmap,
        # legend = :true,
        # phase_time = ts[30],
        # c = :lajolla,
        field_apply = field_apply,
        # clims = (minc, maxc),
        # leg = false,
        title = "",
        frame = :none
    )
    plot!(bearing)
    # plot!(xlim = [-bearing.inner_radius * 0.4,bearing.inner_radius * 0.4], ylims = [-bearing.outer_radius*0.8, -bearing.inner_radius * 0.85])

## Calculate results for all frequencies
    # inds = 1:7
    # ωs = ωs[inds];

    results = map(ωs) do ω

        bd1_for = BoundaryData(ω, bearing, loading_profile; σ = σ)
        basis_order = maximum(abs.(bd1_for.modes)) + 5;
        θs = LinRange(0.0, 2pi, 2basis_order+2)[1:end-1]
    
        bd2_for = BoundaryData(bc2_forward, θs=θs, fields = hcat(0.0.*θs .+ 0.0im, 0.0.*θs .+ 0.0im))

        # bd2_for = BoundaryData(bc2_forward, modes = bd1_for.modes, coefficients =  0.0 .* bd1_for.coefficients)

        sim = BearingSimulation(ω, bearing, bd1_for, bd2_for; 
            method = modal_method
        );
        wave = ElasticWave(sim)

        result_traction = field(wave, bearing, TractionType(); res = 150);
        fs = [f[1] for f in result_traction.field];
        result_pressure = FrequencySimulationResult(fs, result_traction.x, [ω]);
    end;

    # the norm should decrese at some point with an increase in frequency
    [maximum(norm.(r.field)) for r in results]

    inds = eachindex(results);
    # inds = 1:3

    fields = [[f[1] for f in field(r)] for r in results[inds]];
    all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs[inds]);

    gr(size = (400,400))
    i = 3;
    plot(all_results, ωs[i];
        seriestype=:heatmap
        , field_apply = f -> abs(f[1])
    )

## Plot a time video
    ts = ω_to_t(ωs[inds])
    
    ts = LinRange(0,total_time / Z,60)[1:end-1]
    # ts = LinRange(0,ts[end],43)
    time_result = frequency_to_time(all_results; t_vec = ts)

    # maxc = 0.16 .* maximum(norm.(field(time_result)))
    maxc = 0.6 .* maximum(norm.(field(time_result)))
    # maxc = 0.5 .* maximum(norm.(field(time_result)))
    minc = - maxc

    t = ts[2]
    t = ts[end]

    # pyplot()

    anim = @animate for t in ts[1:end]
        plot(time_result, t,
            seriestype=:heatmap,
            # field_apply = f -> real(f[2]),
            clim = (minc, maxc),
            leg = false,
        )
        # scatter!([bearing.inner_radius * cos(θo)], [bearing.inner_radius * sin(θo)])
        plot!(bearing, t)
        plot!(frame = :none, title="", xguide ="",yguide ="")
    end
    
    gif(anim,"docs/images/bearings-time.gif", fps = 5)
    # gif(anim,"docs/images/bearings-time-high-freq.gif", fps = 7)
    # gif(anim,"docs/images/bearings-time-crack.gif", fps = 4)
    # gif(anim,"docs/images/bearings-time-slow.gif", fps = 4)