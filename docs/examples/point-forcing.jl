ω = 4.0
steel = Elasticity(2; ρ = 7.0, cp = 5.0 - 1.1im, cs = 3.5 - 0.6im)
bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius=2.0)

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
fp = θs .* 0 + θs .* 0im |> collect;

# fp[1:3] = 1e5 .* ones(3) - 1e5im .* ones(3)
fp[1] = 1e5

θ0 = θs[1]
x0 = radial_to_cartesian_coordinates([bearing.inner_radius, θ0])

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

## The long way to calculate result = field(wave.pressure, bearing; res = 100)
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

    fs = [field(wave.pressure,x) for x in xs];
    field_mat[inds,:] = fs

    result2 = FrequencySimulationResult(field_mat, x_vec, [ω])

## plot the field in the bearing

using Plots
plot(result,ω;
    seriestype = :contour,
    # seriestype = :heatmap,
    # legend = :false,
    field_apply = real
)
scatter!([x0[1]],[x0[2]])
# savefig("docs/images/bearing-stress-patterns.png")

maxc = maximum(real.(field(result)))
minc = minimum(real.(field(result)))

maxc = min(abs(minc),abs(maxc))
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

maxω = 4ω

ωs = LinRange(0.01,maxω,100)

basis_order = estimate_basisorder(maxω, bearing; tol =1e-5)
basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

# 0.0 and 2pi are the same point
θs = LinRange(0.0,2pi, basis_length + 1)[1:end-1]

# let's create a focused pressure on the inner boundary
fs = θs .* 0 + θs .* 0im |> collect;
fp = θs .* 0 + θs .* 0im |> collect;
fp[1] = 1.0 + 0.0im

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


results = map(eachindex(ωs)) do i
    sim = BearingSimulation(ωs[i], bearing, bd1, bd2)
    wave = ElasticWave(sim)

    field(wave.pressure, bearing; res = 120)
end

pyplot()


fields = [r.field for r in results];
all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs);

## NOTE the amplitude of the field for ωs[1] is too high. It is dominating the Fourier transform. Need to investigate

i = 1;
plot(all_results, ωs[i]; seriestype=:contour)

maxc = 1.5mean(real.(field(all_results)))
# maxc = maxc /5.0
minc = - maxc

anim = @animate for ω in ωs
    plot(all_results, ω,
      seriestype=:contour,
      # seriestype=:heatmap,
      # clim = (minc, maxc),
      leg = false,
    )
    plot!(frame = :none, title="", xguide ="",yguide ="")
end

gif(anim,"docs/images/all-frequencies-point-pressure.gif", fps = 3)

ts = ω_to_t(ωs)
impulse = GaussianImpulse(maxω)
impulse.in_time.(ts)
frequency_impulse = impulse.in_freq.(ωs)

time_result = frequency_to_time(all_results; t_vec = ts, impulse = impulse)

maxc = maximum(abs.(field(time_result)))
# maxc = mean(abs.(field(time_result)))
# maxc = maxc /5.0
minc = - maxc

t = ts[40]
t = ts[140]
t = ts[10]
anim = @animate for t in ts
    plot(time_result, t,
      seriestype=:contour,
      # seriestype=:heatmap,
      clim = (minc, maxc),
      # leg = false,
    )
    plot!(frame = :none, title="", xguide ="",yguide ="")
end

gif(anim,"docs/images/time-point-pressure.gif", fps = 7)
