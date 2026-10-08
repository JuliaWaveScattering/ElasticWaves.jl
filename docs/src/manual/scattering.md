# Scattering from spheres

ElasticWaves.jl extends the package [MultipleScattering.jl](https://github.com/JuliaWaveScattering/MultipleScattering.jl) to elastic waves in three spatial dimensions. This means that all the tools from MultipleScattering.jl, such as `Particle`, `FrequencySimulation`, `run`, and the plot recipes, can be used with an [`Elastic`](@ref) medium. See the [documentation of MultipleScattering.jl](https://juliawavescattering.github.io/MultipleScattering.jl/dev/) for more on these tools.

## Two spheres

First we define the background medium, and the medium of the particles:

```@example scattering
using ElasticWaves, MultipleScattering, Plots

medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 0.6)
inclusion = Elastic(3; ρ = 3.0, cp = 2.0, cs = 1.2)
nothing # hide
```

Next we place two spheres, with radius 1, side by side:

```@example scattering
particles = [
    Particle(inclusion, Sphere([-2.0, 0.0, 0.0], 1.0)),
    Particle(inclusion, Sphere([ 2.0, 0.0, 0.0], 1.0))
]
nothing # hide
```

The source is a plane shear wave which travels along the z-axis and is polarised along the x-axis, see [`plane_z_shear_source`](@ref).

```@example scattering
source = plane_z_shear_source(medium)
sim = FrequencySimulation(particles, source)
nothing # hide
```

To get results we `run` the simulation for an angular frequency `ω` and a region. A region in three dimensions leads to a grid of points on a plane with a constant `y`. Below this is the plane `y = 0`, which passes through the centres of both spheres.

```@example scattering
ω = 1.2
region = Box([[-5.0, 0.0, -6.0], [5.0, 0.0, 6.0]])

result = run(sim, region, ω; res = 121)
nothing # hide
```

The field of the `result` is the displacement, with its three Cartesian components, at each point of the grid. This includes the points inside the spheres.

```@example scattering
field(result)[1]
```

The incident wave has only a displacement in the x direction. The displacement in the z direction is due only to the waves scattered from the spheres, which are both pressure and shear waves.

```@example scattering
circles = [Circle(origin(p)[[1, 3]], outer_radius(p)) for p in particles]

p1 = plot(result, ω; seriestype = :heatmap,
    field_apply = u -> real(u[1]), clims = (-1.8, 1.8)
)
plot!(p1, circles[1], linecolor = :black)
plot!(p1, circles[2], linecolor = :black)
plot!(p1, title = "x-displacement", yguide = "z")

p2 = plot(result, ω; seriestype = :heatmap,
    field_apply = u -> real(u[3]), clims = (-0.9, 0.9)
)
plot!(p2, circles[1], linecolor = :black)
plot!(p2, circles[2], linecolor = :black)
plot!(p2, title = "z-displacement", yguide = "z")

plot(p1, p2, layout = (1, 2), size = (900, 450), margin = 4Plots.mm)
savefig("sphere-scattering.png"); nothing # hide
```
![The displacement around two spheres](sphere-scattering.png)

!!! warning "The axis through the centre of each sphere"
    The vector spherical harmonics used for the displacement are not defined on the line which passes through the centre of a sphere and is parallel to the z-axis. Evaluating the field on one of these lines returns `NaN`. Above we chose `res = 121` so that no point of the grid has `x = -2` or `x = 2`.

## The T-matrix

The scattering from each particle is described by its T-matrix, which relates the coefficients of the wave that excites the particle to the coefficients of the wave scattered from the particle. For an elastic sphere in an elastic medium use [`t_matrix`](@ref):

```@example scattering
basis_order = 5
T = t_matrix(particles[1], medium, ω, basis_order)
size(T)
```

For each spherical harmonic there are three coefficients: one for the pressure potential and two for the shear (Debye) potentials. So the number of coefficients for one particle is

```@example scattering
basisorder_to_basislength(typeof(medium), basis_order)
```

The functions that these coefficients multiply are given by [`regular_basis_function`](@ref) and [`outgoing_basis_function`](@ref). The derivation of the T-matrix is in the notes [ElasticSphereScattering.pdf](https://github.com/JuliaWaveScattering/ElasticWaves.jl/blob/main/docs/theory/ElasticSphereScattering.pdf).
