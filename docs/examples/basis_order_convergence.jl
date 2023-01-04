ω = 0.1
steel = Elasticity(2; ρ = 7.0, cp = 5.0 - 1.1im, cs = 3.5 - 0.6im)
bearing = RollerBearing(medium = steel, inner_radius = 1.0, outer_radius = 2.0)

# this non-dimensional number determines what basis_order is neeeded
kpa = bearing.outer_radius * ω / steel.cp
ksa = bearing.outer_radius * ω / steel.cs

# estimate the largest basis_order that a wave scattered from the inner boundary can be measured at the outer boundary

basis_orders = 4:5:40 |> collect
basis_lengths = basisorder_to_basislength.(Acoustic{Float64,2}, basis_orders)

# 0.0 and 2pi are the same point
θs = LinRange(0.0,2pi, basis_lengths[end] + 1)[1:end-1]

# let's create a focused pressure on the inner boundary
fs = θs .* 0 + θs .* 0im |> collect;
fp = exp.(-6.0 .* (θs .- pi).^2) + θs .* 0im |> collect;

x0 = bearing.inner_radius .* [-1.0,0.0];

bd1 = BoundaryData(TractionBoundary(inner = true); θs = θs, fields = hcat(fp,fs))
bd2 = BoundaryData(TractionBoundary(outer = true); θs = θs, fields = hcat(fs,fs))


sims = [
    BearingSimulation(ω, bearing, bd1, bd2; basis_order = m)
for m in basis_orders]

waves = ElasticWave.(sims);

results = [field(w.pressure, bearing; res = 100) for w in waves];

using Plots
gr(size = (1000,600))
ps = [
    plot(r, ω;
        seriestype=:heatmap
    )
for r in results]

plot(ps...)
savefig("docs/images/basis_convergence.pdf")
