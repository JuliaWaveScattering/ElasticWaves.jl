ω = 4.0
ω = 20000.0
steel = Elastic(2; ρ = 7.0, cp = 5.0 - 0.1im, cs = 3.5 - 0.2im)
steel = Elastic(2; ρ = 7800.0, cp = 5000.0 -0.5im, cs = 3500.0 -0.5im)

bearing = RollerBearing(medium = steel, inner_radius = 1.0, outer_radius = 2.0)

kpa = bearing.outer_radius * ω / steel.cp
bearing.inner_radius * ω / steel.cp
ksa = bearing.outer_radius * ω / steel.cs

# estimate the largest basis_order that a wave scattered from the inner boundary can be measured at the outer boundary

basis_order = 10
# basis_order = estimate_basisorder(ω,bearing; tol =1e-5)
basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

# 0.0 and 2pi are the same point
θs = LinRange(0.0,2pi, basis_length + 1)[1:end-1]

# let's create a focused pressure on the inner boundary
fs = θs .* 0 + θs .* 0im |> collect;
fp = exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im |> collect;

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
inner_field = fouriermodes_to_fields(θs2,sim.boundarydata1.coefficients)

using Plots
plot(θs2,real.(inner_field))
plot!(θs, real.(fp))

wave = ElasticWave(sim);
result = field(wave.potentials[1], bearing; res = 120)

## The long way to calculate result = field(bearing, wave.potentials[1]; res = 100)
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

    fs = [field(wave.potentials[1], x) for x in xs];
    field_mat[inds,:] = fs

    result2 = FrequencySimulationResult(field_mat, x_vec, [ω])

## plot the field in the bearing

ts = (0.0:0.1:2pi) ./ ω

using Plots
gr(size = (400,400))

sim = BearingSimulation(ω, bearing, bd1, bd2)
wave = ElasticWave(sim);
result = field(wave.potentials[1], bearing; res = 420)
# result = field(wave.potentials[2], bearing; res = 420)

field_apply = abs
maxc = 0.8 .* maximum(field_apply.(field(result)))
minc = min(0.0,field_apply(- maxc))

plot(result,ω;
    # seriestype = :contour,
    seriestype = :heatmap,
    # legend = :false,
    phase_time = ts[30],
    # c = :lajolla,
    field_apply = field_apply,
    clims = (minc, maxc),
    leg = false
)
plot!(frame = :none, title="Pressure wave", xguide ="",yguide ="")
# plot!(frame = :none, title="Shear wave", xguide ="",yguide ="")
plot!(Circle(bearing.inner_radius))
plot!(Circle(bearing.outer_radius))
plot!(Circle([r-bearing.inner_radius, 0.0], r))
pp = plot!(Circle(bearing.inner_radius -2r))
# ps = plot!(Circle(bearing.inner_radius -2r))
# scatter!([x0[1]],[x0[2]], lab = "")
# savefig("docs/images/bearing-shear-stress-patterns.png")
# savefig("docs/images/bearing-shear-stress-patterns.pdf")
savefig("docs/images/bearing-pressure-stress-patterns.png")
savefig("docs/images/bearing-pressure-stress-patterns.pdf")

sim = BearingSimulation(ω, bearing, bd1, bd2)
wave = ElasticWave(sim);
result = field(wave.potentials[1], bearing; res = 120)

field_apply = abs
maxc = 0.8 .* maximum(field_apply.(field(result)))
minc = min(0.0,field_apply(- maxc))

plot(result,ω;
    # seriestype = :contour,
    seriestype = :heatmap,
    # legend = :false,
    phase_time = ts[30],
    # c = :lajolla,
    field_apply = field_apply,
    clims = (minc, maxc),
    leg = false
)
plot!(frame = :none, title="Roller bearing pressure", xguide ="",yguide ="")
plot!(Circle(bearing.inner_radius))
plot!(Circle(bearing.outer_radius))
plot!(Circle([r-bearing.inner_radius, 0.0], r))
plot!(Circle(bearing.inner_radius -2r))
# scatter!([x0[1]],[x0[2]], lab = "")
savefig("docs/images/bearing-stress-patterns.png")

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
