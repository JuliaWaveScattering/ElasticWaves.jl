using ElasticWaves

using Test, Statistics, LinearAlgebra, MultipleScattering
using Plots 

# the higher the frequency, the worse the result. This is already a high frequency.
medium = Elastic(2; ρ = 2.0, cp = 2.0 - 0.01im, cs = 1.0 - 0.02im)

Ω = 0.8 # the angular speed should be much smaller than the wavesppeds, otherwise the figure does not make sense

bearing = RollerBearing(medium = medium, 
    inner_radius = 1.2, outer_radius = 2.0, 
    angular_speed = Ω,  
    rollers_inside = true
)

plot(bearing, 0.0)


# Calculate the frequencies to be used for the time video

    # we want enough resolution in time to see the rollers going around. 
    total_time = 2pi / bearing.angular_speed

    # to get the number of frames per cycle:
    frames = 200
    dt = total_time / frames

    # from dt we can calculate what's the maximum frequency we need, and then the range of natural frequencies we should use
    maxω = 2pi / dt
    frequency_order = Int(round(maxω / (bearing.number_of_rollers * bearing.angular_speed)))

    ωms = natural_frequencies(bearing, frequency_order) |> collect

    # need to exclude the zero frequency as is ill posed
    ωs = ωms[2:end]


## Do one example
    ω = ωms[end]

    dr = bearing.outer_radius - bearing.inner_radius
    kp = (ω / medium.cp)
    kp * dr

    loading_basis_order = 15;    
    loading_θs = LinRange(0.0, 2pi, 2loading_basis_order+2)[1:end-1]

    # loading profile
    θo = 3pi/2;
    fp_loading = 0.3 .- exp.(-1.0 .* (loading_θs .- θo).^2) + loading_θs .* 0im; 
    # fp_loading = -0.3 .+ exp.(-1.0 .* (loading_θs .- θo).^2) + loading_θs .* 0im; 
    fs_loading = 0.4 .* fp_loading;


    plot(loading_θs, real.(fp_loading))

    bc1_forward = TractionBoundary(inner=true)
    bc2_forward = TractionBoundary(outer=true)

    loading_profile = BoundaryData(bc1_forward, θs=loading_θs, fields = hcat(fp_loading,fs_loading))

    basis_order = 15;
    θs = LinRange(0.0, 2pi, 2basis_order+2)[1:end-1]

    bd1_for = BoundaryData(ω, bearing, loading_profile)
    bd2_for = BoundaryData(bc2_forward, θs=θs, fourier_modes = 0.0 .* bd1_for.fourier_modes)

    modal_method = ModalMethod(tol = 1e-9, only_stable_modes = true)
    forward_sim = BearingSimulation(ω, bearing, bd1_for, bd2_for; 
        method = modal_method,
        nondimensionalise = true);

    wave = ElasticWave(forward_sim);

    result = field(wave.potentials[1], bearing; res = 100)

    field_apply = real
    maxc = 0.5 .* maximum(field_apply.(field(result)))
    minc = min(0.0,field_apply(- maxc))

    plot(result,ω;
        seriestype = :heatmap,
        # legend = :true,
        # phase_time = ts[30],
        # c = :lajolla,
        field_apply = field_apply,
        clims = (minc, maxc),
        leg = false
    )
    plot!(bearing)

## Calculate results for all frequencies

    results = map(ωs) do ω

        bd1_for = BoundaryData(ω, bearing, loading_profile)

        basis_order = basislength_to_basisorder(PhysicalMedium{2,1},size(bd1_for.fourier_modes,1))
        
        # θs = LinRange(0.0, 2pi, 2basis_order+2)[1:end-1]

        bd2_for = BoundaryData(bc2_forward, fourier_modes = 0.0 .* bd1_for.fourier_modes)

        sim = BearingSimulation(ω, bearing, bd1_for, bd2_for; 
            method = modal_method,
            nondimensionalise = true
        );
        wave = ElasticWave(sim)

        # res = field(wave, bearing, TractionType(); res = 70)

        # # scale the potential to match the units of stress
        scale = medium.ρ * ω^2

        coes = wave.potentials[1].coefficients;
        # coes[2,:] .= 0.0 + 0.0im

        potential = HelmholtzPotential{2}(wave.potentials[1].wavespeed, wave.potentials[1].wavenumber, scale .* coes)

        res = field(potential, bearing; res = 100)
    end

    fields = [[f[1] for f in field(r)] for r in results];
    all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs);

    i = 3;
    plot(all_results, ωs[i];
        seriestype=:heatmap
        , field_apply = f -> real(f[1])
    )

## Plot a time video
    ts = ω_to_t(ωs)
    time_result = frequency_to_time(all_results; t_vec = ts)

    maxc = 0.14 .* maximum(norm.(field(time_result)))
    minc = - maxc

    t = ts[4]

    # pyplot()

    anim = @animate for t in ts[1:end]
        plot(time_result, t,
            seriestype=:heatmap,
            # field_apply = f -> real(f[2]),
            clim = (minc, maxc),
            leg = false,
        )
        plot!(bearing, t)
        plot!(frame = :none, title="", xguide ="",yguide ="")
    end
    
    gif(anim,"docs/images/bearings-time.gif", fps = 8)