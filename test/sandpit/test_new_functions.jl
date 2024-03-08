using ElasticWaves
using MultipleScattering
using LinearAlgebra
using Plots


basis_order = 11;
numberofsensors = 2
basis_length = 2*basis_order + 1

#Friction coefficient
#μ=1
θs = LinRange(0, 2pi, 20*basis_length + 2)[1:end-1] |>collect
#θs = LinRange(0, 2pi, 401)[1:end-1] 
#θ2s = LinRange(0, 2pi, 4*basis_length + 1)[1:end-1]
θ2s = LinRange(0, 2pi, 2000)[1:end-1] 
θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]


#Friction coefficient


#Properties of the bearing

steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0,friction_coefficient=0.5)
bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0, number_of_rollers=11, angular_speed=1000.0)
μ=bearing.medium.friction_coefficient
Ω=bearing.angular_speed

Z=bearing.number_of_rollers

i=10
n_order=basis_order
ωs=[n*Z*Ω for n in 0:n_order] 


#Angular velocity of the bearings


bc1_forward = TractionBoundary(inner=true)
bc2_forward = TractionBoundary(outer=true)

bc1_inverse = DisplacementBoundary(outer=true)
bc2_inverse = TractionBoundary(outer=true)

bd1_forward=point_contact_boundary_data(θs,bearing,bc1_forward,basis_order)


f0=0*θs

bd2_forward=bd2_forward = BoundaryData(bc2_forward,θs=θs, fields=hcat(f0,f0))

modal_method=ModalMethod(tol=1e-6,only_stable_modes=true)

sim = BearingSimulation(ωs[i], bearing, bd1_forward, bd2_forward
; method=modal_method,nondimensionalise=true)



wave = ElasticWave(sim)



#x_outer are the coordinates where our sensors are
x_outer=[radial_to_cartesian_coordinates([bearing.outer_radius,θ]) for θ in θs_inv ]

#x2_inner are the cooordinates that we will calculate the forces that the fields of our solutions generates
x2_inner = [
    radial_to_cartesian_coordinates([bearing.inner_radius, θ])
for θ in θ2s]


x2_outer = [
    radial_to_cartesian_coordinates([bearing.outer_radius, θ])
for θ in θ2s]



traction_for= [traction(wave,x) for x in x2_inner]
traction_for = hcat(traction_for...) |> transpose |> collect

plot(θ2s,real.(traction_for[:,1]))

q0=1
ϵ=0.4


load = q0*(1 .-(1/(2*ϵ)).*(1 .- cos.(θs)) + θs .* 0im).^(10) 

load=load.^(1/9) |> real

plot(θs,real.(load))


frequency_order=i

load=LoadingProfile(bc1_forward, θs=θs, fields=hcat(load, μ*load))

bd1_forward=BoundaryData(load, bearing, basis_order ,frequency_order) 


bd2_forward=bd2_forward = BoundaryData(bc2_forward,θs=θs, fields=hcat(f0,f0))

modal_method=ModalMethod(tol=1e-6,only_stable_modes=true)

sim = BearingSimulation(ωs[i], bearing, bd1_forward, bd2_forward
; method=modal_method,nondimensionalise=true)



wave = ElasticWave(sim)



#x_outer are the coordinates where our sensors are
x_outer=[radial_to_cartesian_coordinates([bearing.outer_radius,θ]) for θ in θs_inv ]

#x2_inner are the cooordinates that we will calculate the forces that the fields of our solutions generates
x2_inner = [
    radial_to_cartesian_coordinates([bearing.inner_radius, θ])
for θ in θ2s]


x2_outer = [
    radial_to_cartesian_coordinates([bearing.outer_radius, θ])
for θ in θ2s]



traction_for= [traction(wave,x) for x in x2_inner]
traction_for = hcat(traction_for...) |> transpose |> collect

plot(θ2s,real.(traction_for[:,1]))
plot!(θs,real.(load))
