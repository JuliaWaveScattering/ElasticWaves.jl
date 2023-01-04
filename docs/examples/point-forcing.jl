## Simulation in time

steel = Elasticity(2; ρ = 7.0, cp = 5.0 - 1.1im, cs = 3.5 - 0.6im)
bearing = RollerBearing(medium = steel, inner_radius = 1.0, outer_radius = 2.0)

# time for wave to go from inner boundary to outer boundary
Δt = (bearing.outer_radius - bearing.inner_radius) / steel.cp |> real

# need at as many samples as possible within this time frame
dt = Δt / 4

# from dt we can calculate what's the maximum frequency we need
maxω = 2pi / dt
ωs = LinRange(0.3,maxω,500)

kpa = bearing.outer_radius * maxω / steel.cp
ksa = bearing.outer_radius * maxω / steel.cs

kpa = bearing.inner_radius * ωs[1] / steel.cp
ksa = bearing.inner_radius * ωs[1] / steel.cs

basis_order = 12
basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

# 0.0 and 2pi are the same point
θs = LinRange(0.0,2pi, basis_length + 1)[1:end-1]

# let's create a focused pressure on the inner boundary
fs = θs .* 0 + θs .* 0im |> collect;
fp = exp.(-6.0 .* (θs .- pi).^2) + θs .* 0im |> collect;

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
res = field(wave.pressure, bearing; res = 120)

# @time field(wave.pressure, bearing; res = 120);
# @time field(wave, bearing, DisplacementType(); res = 120);

@time field(wave.pressure, rand(2));
@time field(wave, rand(2), TractionType());
@time field(wave, rand(2), DisplacementType());

plot(res,ωs[i]; seriestype=:heatmap)
plot!(Circle(bearing.inner_radius))

field_type = TractionType()
res = field(wave, bearing, field_type; res = 120)

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

    # fieldtype = TractionType()
    # res = field(wave, bearing, fieldtype; res = 120)
    res = field(wave.pressure, bearing; res = 120)
    # plot(res, ωs[i]; seriestype=:contour)
end

# pyplot()

i0 = 3
fields = [r.field for r in results[i0:end]];
all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs[i0:end]);

## NOTE the amplitude of the field for ωs[1] is too high. It is dominating the Fourier transform. Need to investigate

i = 5;
plot(all_results, ωs[i];
    seriestype=:heatmap,
    field_apply = f -> real(f[1])
)
# plot(all_results, ωs[i]; seriestype=:contour)

maxc1 = mean(abs.(real.(field(all_results))))
maxc2 = maximum(abs.(real.(field(all_results))))
maxc = (30 * maxc1 + maxc2) / 31
minc = - maxc

anim = @animate for ω in ωs
    plot(all_results, ω,
      # seriestype=:contour,
      seriestype=:heatmap,
      clim = (minc, maxc),
      leg = false,
    )
    plot!(frame = :none, title="", xguide ="",yguide ="")
end

gif(anim,"docs/images/all-frequencies-point-pressure.gif", fps = 4)

ts = ω_to_t(ωs[i0:end])
impulse = GaussianImpulse(maxω)
plot(ts, impulse.in_time.(ts), xlims = (0.0,0.2))

frequency_impulse = impulse.in_freq.(ωs)
plot(ωs, real.(frequency_impulse))

time_result = frequency_to_time(all_results; t_vec = ts, impulse = impulse)

maxc = maximum(abs.(field(time_result)))
# maxc = mean(abs.(field(time_result)))
# maxc = maxc /5.0
minc = - maxc

t = ts[240]
t = ts[140]
t = ts[20]

# pyplot()
anim = @animate for t in ts
    plot(time_result, t,
      # seriestype=:heatmap,
      seriestype=:contour,
      levels = 24,
      clim = (minc, maxc),
      leg = false,
    )
    plot!(frame = :none, title="", xguide ="",yguide ="")
end

gif(anim,"docs/images/time-point-pressure.gif", fps = 7)
