using MultipleScattering
using Plots
using LinearAlgebra
using SpecialFunctions
using ElasticWaves
using Statistics
using LaTeXStrings
using Accessors
import StaticArrays: SVector

medium = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
r1 = 1.0
r2 = 1.5
bearing = RollerBearing(
    medium = medium, 
    inner_radius = r1, outer_radius = r2
)

# we want wavelengths that go from about 1/10 of the bearing thickness radius to about 2 times the bearing radius
krs = 0.1:0.2:1.0
ks = krs ./ bearing.outer_radius
ωs = real.(ks .* medium.cp)

# Source location
θo = rand(0.0:0.01:2pi)
ro = rand(bearing.inner_radius:0.01:bearing.outer_radius)
xo, yo = radial_to_cartesian_coordinates([ro,θo])
α0=2*rand()
# xo = - bearing.inner_radius / 2 - bearing.outer_radius / 2;
# yo = 0.0;

source_location = [xo,yo]
p_map = SourceMap(source_location |> transpose, [α0])
s_map = SourceMap(source_location |> transpose, [0.0])

basis_order = 4
modes = -basis_order:basis_order |> collect

## The forward problem

# boundarytype = TractionBoundary(inner=true);
displacements_true = map(ωs) do ω 
    # I removed the variable 'θs' from the functions below. Best to just calculate the modes
    bd1 = boundary_data(ω, bearing, TractionBoundary(inner=true), p_map, s_map, modes)
    bd2 = boundary_data(ω, bearing, TractionBoundary(outer=true), p_map, s_map, modes)

    # need to use the minus of the traction for the homogeneous problem, so that the total sum of traction will be zero
    @set bd1.coefficients = - bd1.coefficients
    @set bd2.coefficients = - bd2.coefficients
    method = ModalMethod(modes=modes)

    sim = BearingSimulation(ω, method, bearing, bd1, bd2)
    wave = ElasticWave(sim)

    
    homo_displace = BoundaryData(DisplacementBoundary(outer=true), bearing.outer_radius, wave).coefficients

    source_displace = boundary_data(ω, bearing, DisplacementBoundary(outer=true), p_map, s_map, modes).coefficients

    total_displace = homo_displace + source_displace

    # Add an error of 1%?
    error = 0.01 .* (rand(Float64, size(total_displace)) .- 1.0) .* mean(abs.(total_displace))

    return total_displace + error 
end



## The Inverse problem
# Now we repeat the process above for sources all over the place, and for each of them compare the output displacement
res = 30

inner_circle = Circle(bearing.inner_radius)
outer_circle = Circle(bearing.outer_radius)

x_vec, inds = points_in_shape(outer_circle; res = res, exclude_region = inner_circle)
# x_vec is a square grid of points and x_vec[inds] are the points in the region.

locations = x_vec[inds]

αs=LinRange(0.1,2.0,20) |>collect



fs_α = [zeros(length(locations)) for i in eachindex(αs)]


for i in eachindex(αs)

fs_α[i] = map(locations) do x
    p_map = SourceMap(x |> transpose, [αs[i]])
    s_map = SourceMap(x |> transpose, [0.0])

    displacements = map(ωs) do ω 
        bd1 = boundary_data(ω, bearing, TractionBoundary(inner=true), p_map, s_map, modes) 
        bd2 = boundary_data(ω, bearing, TractionBoundary(outer=true), p_map, s_map, modes)
        @set bd1.coefficients = - bd1.coefficients
        @set bd2.coefficients = - bd2.coefficients

        method = ModalMethod(modes=modes)
    
        sim = BearingSimulation(ω, method, bearing, bd1, bd2)
        wave = ElasticWave(sim)
    
        homo_displace = BoundaryData(DisplacementBoundary(outer=true), bearing.outer_radius, wave).coefficients

        source_displace = boundary_data(ω, bearing, DisplacementBoundary(outer=true), p_map, s_map, modes).coefficients

        return homo_displace + source_displace
    end

    errors = [
    norm(displacements[j] - displacements_true[j]) 
    for j in eachindex(displacements)]

    return 1 ./ norm(errors)
end 

end

fs_α=hcat(fs_α...)

idx=findmax(fs_α)[2][2]

field_mat = [0.0 + 0.0im for x in x_vec]
field_mat[inds] = fs_α[:,idx]

result = FrequencySimulationResult(reshape(field_mat, :, 1), x_vec, [ωs[5]])

field_apply = abs
maxc = 0.8 .* maximum(field_apply.(field(result)))
minc = min(0.0,field_apply(- maxc))

plot(result,ωs[1];
    seriestype = :heatmap,
    c = :inferno,
    field_apply = field_apply,
    clims = (minc, maxc),
)

scatter!([source_location[1]],[source_location[2]], lab = "true source")
plot!(bearing)

α_density = [maximum(fs_α[:,i]) for i in eachindex(αs)] 
plot(αs, α_density)
sum(α_density)
sum(fs_α)