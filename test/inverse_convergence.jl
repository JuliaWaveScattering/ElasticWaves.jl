## We define the inverse problem as taking only information from one boundary, such as the outer boundary, and then predicting fields on the other boundary.
#include("/Users/user/.julia/packages/ElasticWaves/REv3o/test/runtests.jl")

# include("/home/mep21jjk/.julia/packages/ElasticWaves/REv3o/test/runtests.jl")
include("../src/ElasticWaves.jl")

ω = 10.0
# steel = Elasticity(2; ρ = 7.0, cp = 5.0, cs = 3.5)
steel = Elasticity(2; ρ = 7.0, cp = 5.0 - 0.01im, cs = 3.5 - 0.01im)
# steel = Elasticity(2; ρ = 7.0, cp = 5.0, cs = 3.5 )
bearing = RollerBearing(medium = steel, inner_radius = 1.0, outer_radius = 1.1)


kpa = bearing.outer_radius * ω / steel.cp
bearing.inner_radius * ω / steel.cp
ksa = bearing.outer_radius * ω / steel.cs


## Lets define the matrix A for the forward problem 

basis_order = 20
basis_order = 4
ωs = 20.0:0.02:500; 
ω = ωs[2]
data = map(ωs) do ω 
    M = boundarycondition_system(ω, bearing, TractionBoundary(inner=true), TractionBoundary(outer=true), basis_order) 
    N = boundarycondition_system(ω, bearing, DisplacementBoundary(outer=true), TractionBoundary(outer=true), basis_order)
    A = N * inv(M)
    cond(A) 
    # svd(M).S[end] 
end

using Plots 
plot(ωs,data)

svd(M[1:2,1:end]).S
svd(M[3:4,1:end]).S

svd(N).S

# A takes the boundary data from the forward problem and return the boundary data for the inverse problem
A = N * inv(M)
svd(A).S


## Test the complete discrete inverse problem
basis_order = 10
basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

numberofsensors = basis_length - 11;
numberofsensors = 5;

θs = LinRange(0.0, 2pi, basis_length + 1)[1:end-1]
θ2s = LinRange(0.0, 2pi, 4basis_length + 1)[1:end-1]
θs_inv = LinRange(0.0, 2pi, numberofsensors + 1)[1:end-1]

fs = θs .* 0 + θs .* 0im |> collect;
fp = exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im |> collect;

using Plots

bd_traction_inner = BoundaryData(
    TractionBoundary(inner = true);
    θs = θs,
    fields = hcat(fp,fs)
)
bd_traction_outer = BoundaryData(
    TractionBoundary(outer = true);
    θs = θs,
    fields = hcat(fs,fs)
)

# let's have a look at the modes that were calculated during the Bearing. This is the field we will actual approximate

sim = BearingSimulation(ω, bearing, bd_traction_inner, bd_traction_outer)
inner_field = fouriermodes_to_fields(θ2s,sim.boundarydata1.fourier_modes)

plot(θs,real.(fp))
plot!(θ2s,real.(inner_field[:,1]), linestyle = :dash)

wave = ElasticWave(sim)

x_outer = [
    radial_to_cartesian_coordinates([bearing.outer_radius, θ])
for θ in θs_inv]

x_inner = [
    radial_to_cartesian_coordinates([bearing.inner_radius, θ])
for θ in θs_inv]

traction_outer = [traction(wave,x) for x in x_outer];
traction_outer = hcat(traction_outer...) |> transpose |> collect

plot(θs_inv,real.(traction_outer[:,1]))

traction_inner = [traction(wave,x) for x in x_inner];
traction_inner = hcat(traction_inner...) |> transpose |> collect

plot(θs_inv,real.(traction_inner[:,1]))
plot!(θ2s,real.(inner_field[:,1]), linestyle = :dash)

plot(θs_inv,real.(traction_inner[:,2]))

displacement_outer = [displacement(wave,x) for x in x_outer];
displacement_outer = hcat(displacement_outer...) |> transpose |> collect;

plot(θs_inv,real.(displacement_outer[:,1]))

bd_displacement_outer = BoundaryData(
    DisplacementBoundary(outer=true);
    θs = θs_inv,
    fields = displacement_outer
)

inverse_sim = BearingSimulation(ω, bearing, bd_traction_outer, bd_displacement_outer)

outer_field = fouriermodes_to_fields(θ2s, inverse_sim.boundarydata2.fourier_modes)

plot(θs,real.(fp))
plot!(θ2s,real.(inner_field[:,1]), linestyle = :dash)

inverse_wave = ElasticWave(inverse_sim);

# norm(inverse_wave.pressure.coefficients - wave.pressure.coefficients) / norm(wave.pressure.coefficients)
# norm(inverse_wave.shear.coefficients - wave.shear.coefficients) / norm(wave.shear.coefficients)

x2_inner = [
    radial_to_cartesian_coordinates([bearing.inner_radius, θ])
for θ in θ2s]

traction_inv = [traction(inverse_wave,x) for x in x2_inner];
traction_inv = hcat(traction_inv...) |> transpose |> collect

wave_traction = [traction(wave,x) for x in x2_inner];
wave_traction = hcat(wave_traction...) |> transpose |> collect

# plot(θs2,real.(inner_field))
plot(θs,real.(fp))
plot!(θ2s,real.(wave_traction[:,1]), linestyle = :dash)
plot!(θ2s,real.(traction_inv[:,1]),seriestype=:scatter)
plot!(θ2s,real.(traction_inv[:,2]),seriestype=:scatter)
# plot!(θs_inv,real.(fs_pred))

# plot!(θs_inv,real.(fp_pred),seriestype=:scatter)
# plot!(θs_inv,real.(fs_pred))
