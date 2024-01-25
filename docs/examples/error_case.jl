using ElasticWaves

ωs = [10.0,4e4,7e4]

steel = Elastic(2; ρ = 7800.0, cp = 5000.0 -0.5im, cs = 3500.0 -0.5im)
steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0)

# this non-dimensional number determines what basis_order is neeeded
kpas = bearing.outer_radius .* ωs ./ steel.cp
ksas = bearing.inner_radius .* ωs ./ steel.cs

basis_order = 10

basis_length= 2*basis_order+1

forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

bd1 = BoundaryData(TractionBoundary(inner=true); fourier_modes=forcing_modes[:, 1:2])
bd2 = BoundaryData(TractionBoundary(outer=true); fourier_modes=forcing_modes[:, 3:4])

sims = [BearingSimulation(ω, bearing, bd1, bd2) for ω in ωs];
waves = [ElasticWave(s) for s in sims];

θs = -pi:0.1:pi; θs |> length;
exps = [
        exp(im * θ * m)
    for θ in θs, m = -basis_order:basis_order];

    # the given forcing
forcing = exps * forcing_modes;
forcing_inner = forcing[:,1:2];
forcing_outer = forcing[:,3:4];

    # the predicted forcing
xs = [bearing.inner_radius .* [cos(θ),sin(θ)] for θ in θs];
forcing_inners = [
        hcat(([traction(w, x) for x in xs])...) |> transpose
    for w in waves];

xs = [bearing.outer_radius .* [cos(θ),sin(θ)] for θ in θs];
    forcing_outers = [
        hcat(([traction(w, x) for x in xs])...) |> transpose
    for w in waves];

errors = [maximum(abs.(forcing_inner - f)) for f in forcing_inners];
