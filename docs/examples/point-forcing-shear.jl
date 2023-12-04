## Simulation in time

steel = Elastic(2; ρ = 7.0, cp = 5.0 - 1.1im, cs = 3.5 - 0.6im)
steel = Elastic(2; ρ = 7.0, cp = 5.0 - 0.05im, cs = 3.5 - 0.1im)
bearing = RollerBearing(medium = steel, inner_radius = 1.0, outer_radius = 2.0)

# time for wave to go from inner boundary to outer boundary
Δt = (bearing.outer_radius - bearing.inner_radius) / steel.cp |> real

# need at as many samples as possible within this time frame
dt = Δt / 4

# from dt we can calculate what's the maximum frequency we need
maxω = 2pi / dt
ωs = LinRange(0.2,maxω,500)

kpa = bearing.outer_radius * maxω / steel.cp
ksa = bearing.outer_radius * maxω / steel.cs

kpa = bearing.inner_radius * ωs[1] / steel.cp
ksa = bearing.inner_radius * ωs[1] / steel.cs

basis_order = 24
basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

# 0.0 and 2pi are the same point
θs = LinRange(0.0,2pi, basis_length + 1)[1:end-1]

# let's create a focused pressure on the inner boundary
fs = θs .* 0 + θs .* 0im |> collect;
fp = exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im |> collect;

using Plots
plot(θs,real.(fp))

bd1 = BoundaryData(
    TractionBoundary(inner = true);
    fields = hcat(fp,fs),
    θs = θs,
)
bd2 = BoundaryData(
    TractionBoundary(outer = true);
    fields = hcat(fs,fs),
    θs = θs,
)

i = 2
i = 100
# sim = BearingSimulation(ωs[i], bearing, bd1, bd2; tol = 1e-9)
sim = BearingSimulation(ωs[i], bearing, bd1, bd2; basis_order = basis_order, tol = 1e-9)

wave = ElasticWave(sim)
res = field(wave.potentials[2], bearing; res = 120)

plot(res,ωs[i]; seriestype=:heatmap)
plot!(Circle(bearing.inner_radius))

field_type = TractionType()
res = field(wave, bearing, field_type; res = 80)

# plot the radial traction
plot(res,ωs[i]; seriestype=:heatmap, field_apply = f -> real(f[1]))
plot!(Circle(bearing.inner_radius))
plot!(Circle(bearing.outer_radius))

# plot the theta traction
plot(res,ωs[i]; seriestype=:heatmap, field_apply = f -> real(f[2]))
plot!(Circle(bearing.inner_radius))
plot!(Circle(bearing.outer_radius))

θ2s = LinRange(0.0,2pi, 3basis_length + 1)[1:end-1]

xs = [bearing.inner_radius .* [cos(θ), sin(θ)] for θ in θ2s];
trs = [traction(wave, x)[1] for x in xs]

plot(θ2s,real.(trs))
plot!(θs,real.(fp), linestyle = :dash)

results = map(eachindex(ωs)) do i
    sim = BearingSimulation(ωs[i], bearing, bd1, bd2;
        basis_order = basis_order)
    wave = ElasticWave(sim)

    # res = field(wave, bearing, TractionType(); res = 70)

    # scale the potential to match the units of stress
    scale = steel.ρ * ωs[i]^2

    potential = HelmholtzPotential{2}(wave.potentials[2].wavespeed, wave.potentials[2].wavenumber, scale .* wave.potentials[2].coefficients)

    res = field(potential, bearing; res = 120)
    # plot(res, ωs[i]; seriestype=:contour)
end

i0 = 1
fields = [[f[1] for f in field(r)] for r in results[i0:end]];
all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs[i0:end]);

## NOTE the amplitude of the field for ωs[1] is too high. It is dominating the Fourier transform. Need to investigate

i = 1;
plot(all_results, ωs[i];
    seriestype=:heatmap
    # , field_apply = f -> real(f[1])
)
# plot(all_results, ωs[i]; seriestype=:contour)

# maxc1 = mean(norm.(real.(field(all_results))))
# maxc2 = maximum(norm.(real.(field(all_results))))
# maxc = (30 * maxc1 + maxc2) / 31
# maxc = (2 * maxc1 + maxc2) / 3
# minc = - maxc
#
# anim = @animate for ω in ωs
#     plot(all_results, ω,
#       # seriestype=:contour,
#       seriestype=:heatmap,
#       # field_apply = f -> real(f[2]),
#       clim = (minc, maxc),
#       leg = false,
#     )
#     plot!(frame = :none, title="", xguide ="",yguide ="")
# end
#
# gif(anim,"docs/images/all-frequencies-point-pressure.gif", fps = 4)

ts = ω_to_t(ωs[i0:end])
impulse = GaussianImpulse(maxω)
plot(ts, impulse.in_time.(ts), xlims = (0.0,0.2))

frequency_impulse = impulse.in_freq.(ωs)
plot(ωs, real.(frequency_impulse))

time_result = frequency_to_time(all_results; t_vec = ts, impulse = impulse)

maxc = 0.33 .* maximum(norm.(field(time_result)))
# maxc = mean(abs.(field(time_result)))
# maxc = maxc /5.0
minc = - maxc

t = ts[240]
t = ts[20]
t = ts[50]
t = ts[1]

# pyplot()
r = bearing.inner_radius / 4

anim = @animate for t in ts[1:80]
    plot(time_result, t,
      seriestype=:heatmap,
      # seriestype=:contour,
      # field_apply = f -> real(f[2]),
      # levels = 30,
      clim = (minc, maxc),
      leg = false,
    )
    plot!(frame = :none, title="", xguide ="",yguide ="")
    plot!(Circle(bearing.inner_radius))
    plot!(Circle(bearing.outer_radius))
    plot!(Circle([r-bearing.inner_radius, 0.0], r))
    plot!(Circle(bearing.inner_radius -2r))
end

gif(anim,"docs/images/time-point-shear.gif", fps = 7)
