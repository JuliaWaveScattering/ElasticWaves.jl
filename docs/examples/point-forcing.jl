ω = 4.0
steel = Elasticity(2; ρ = 7.0, cp = 5.0 - 1.1im, cs = 3.5 - 0.6im)
bearing = RollerBearing(medium = steel, inner_radius = 1.0, outer_radius = 2.0)

# this non-dimensional number determines what basis_order is neeeded
kpa = bearing.outer_radius * ω / steel.cp
bearing.inner_radius * ω / steel.cp
ksa = bearing.outer_radius * ω / steel.cs

# estimate the largest basis_order that a wave scattered from the inner boundary can be measured at the outer boundary

basis_order = estimate_basisorder(ω,bearing; tol =1e-5)
basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

# 0.0 and 2pi are the same point
θs = LinRange(0.0,2pi, basis_length + 1)[1:end-1]

# let's create a focused pressure on the inner boundary
fs = θs .* 0 + θs .* 0im |> collect;
fp = exp.(-6.0 .* (θs .- pi).^2) + θs .* 0im |> collect;

# fp[1:3] = 1e5 .* ones(3) - 1e5im .* ones(3)
# fp[1] = 1e5

# θ0 = θs[1]
# x0 = radial_to_cartesian_coordinates([bearing.inner_radius, θ0])
x0 = bearing.inner_radius .* [-1.0,0.0];

bd1 = BoundaryData(TractionBoundary(inner = true); θs = θs, fields = hcat(fp,fs))
bd2 = BoundaryData(TractionBoundary(outer = true); θs = θs, fields = hcat(fs,fs))

sim = BearingSimulation(ω, bearing, bd1, bd2)

# let's have a look at the modes that were calculated during the Bearing. This is the field we will actual approximate
θs2 = LinRange(0.0,2pi,100)
inner_field = fouriermodes_to_fields(θs2,sim.boundarydata1.fourier_modes)

using Plots
plot(θs2,real.(inner_field))
plot!(θs, real.(fp))

wave = ElasticWave(sim);
result = field(wave.pressure, bearing; res = 120)

## The long way to calculate result = field(bearing, wave.pressure; res = 100)
    inner_circle = Circle(bearing.inner_radius)
    outer_circle = Circle(bearing.outer_radius)

    res = 120
    x_vec, inds = points_in_shape(outer_circle; res = res,
        exclude_region = inner_circle
    )
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.

    scatter([x[1] for x in x_vec], [x[2] for x in x_vec])
    scatter!([x[1] for x in x_vec[inds]], [x[2] for x in x_vec[inds]])

    xs = x_vec[inds]
    field_mat = zeros(Complex{Float64},length(x_vec), 1) # change 1 to number of different frequencies

    fs = [field(wave.pressure, x) for x in xs];
    field_mat[inds,:] = fs

    result2 = FrequencySimulationResult(field_mat, x_vec, [ω])

## plot the field in the bearing

ts = (0.0:0.1:2pi) ./ ω

using Plots
plot(result,ω;
    seriestype = :contour,
    # seriestype = :heatmap,
    # legend = :false,
    phase_time = ts[30],
    field_apply = real
)
scatter!([x0[1]],[x0[2]])
# savefig("docs/images/bearing-stress-patterns.png")

maxc = maximum(real.(field(result)))
minc = minimum(real.(field(result)))

maxc = mean([abs(minc),abs(maxc)])
minc = - maxc

anim = @animate for t in ts
    plot(result,ω;
        seriestype = :contour,
        # seriestype = :heatmap,
        phase_time=t, clim=(minc,maxc),
        c=:balance
    )
    plot!(colorbar=false, title="",axis=false, xguide ="",yguide ="")
end

gif(anim,"docs/images/bearing-point-pressure-2.gif", fps = 7)
# gif(anim,"../images/bearing-point-pressure.gif", fps = 7)


## Simulation in time

steel = Elasticity(2; ρ = 7.0, cp = 5.0 - 1.1im, cs = 3.5 - 0.6im)
bearing = RollerBearing(medium = steel, inner_radius = 1.0, outer_radius = 2.0)

# time for wave to go from inner boundary to outer boundary
Δt = (bearing.outer_radius - bearing.inner_radius) / steel.cp |> real

# need at as many samples as possible within this time frame
dt = Δt / 4

# from dt we can calculate what's the maximum frequency we need
maxω = 2pi / dt
ωs = LinRange(0.4,maxω,500)

ω_to_t(ωs)

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

plot(res,ωs[i]; seriestype=:heatmap)
plot!(Circle(bearing.inner_radius))

fieldtype = TractionType()
res = field(wave, bearing, fieldtype; res = 120)

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

    fieldtype = TractionType()
    res = field(wave, bearing, fieldtype; res = 120)
    # res = field(wave.pressure, bearing; res = 120)
    # plot(res, ωs[i]; seriestype=:contour)
end

# pyplot()

i0 = 1
fields = [r.field for r in results[i0:end]];
all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs[i0:end]);

## NOTE the amplitude of the field for ωs[1] is too high. It is dominating the Fourier transform. Need to investigate

i = 1;
plot(all_results, ωs[i];
    seriestype=:heatmap,
    field_apply = f -> real(f[1])
)
# plot(all_results, ωs[i]; seriestype=:contour)

maxc = 2mean(abs.(real.(field(all_results))))
# maxc = maxc /5.0
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

gif(anim,"docs/images/all-frequencies-point-pressure-2.gif", fps = 4)

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
