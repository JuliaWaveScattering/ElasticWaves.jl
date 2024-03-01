using ElasticWaves
using MultipleScattering
using Statistics
using Plots
using FFTW
using LinearAlgebra

#This code do the same of the Zrollers_contact but we extract only one frequency.

basis_order = 5;
numberofsensors = 11
basis_length = 2*basis_order + 1

#Friction coefficient

μ=0.5
#μ=1
θs = LinRange(0, 2pi, 20*basis_length + 1)[1:end-1]
#θs = LinRange(0, 2pi, 401)[1:end-1] 
#θ2s = LinRange(0, 2pi, 4*basis_length + 1)[1:end-1]
θ2s = LinRange(0, 2pi, 2000)[1:end-1] 
θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]


#Friction coefficient


#Properties of the bearing

steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0, number_of_rollers=10.0)

Z=bearing.number_of_rollers |>Int

#Angular velocity of the bearings

Ω=10

m=10
#frequencies
n_order=basis_order+ Z*m |>Int



ω=m*Z*Ω |>Float64

size=2*n_order+1 |>Int

fourier_coef_p=zeros(ComplexF64,size)

q0=1
ϵ=0.4

load = q0*(1 .-(1/(2*ϵ)).*(1 .- cos.(θs)) + θs .* 0im).^(10) 

load=load.^(1/9)

plot(θs,real.(load))

fouter=load*0

bd1 = BoundaryData(TractionBoundary(inner=true), θs=θs, fields = hcat(load,fouter))
bd1 = fields_to_fouriermodes(bd1, basis_order)

modes=bd1.fourier_modes[:,1]

fourier_coef_s=0*fourier_coef_p

f0=0*fourier_coef_p

for n in -basis_order:basis_order
#    cndash=modes[basis_order+1 + n]
    fourier_coef_p[Int(m*Z) + basis_order + 1 + n] = Z*modes[basis_order+1 + n]/2pi 
end

bd1_forward =  BoundaryData(bc1_forward,θs=θs, fourier_modes=hcat(fourier_coef_p,fourier_coef_s))
    
bd2_forward = BoundaryData(bc2_forward,θs=θs, fourier_modes=hcat(f0,f0))







sim = BearingSimulation(ω, bearing, bd1_forward, bd2_forward; basis_order = n_order)


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




traction_for_ω= [traction(wave,x) for x in x2_inner]
traction_for_ω = hcat(traction_for_ω...) |> transpose |> collect

plot(θ2s,real.(traction_for_ω[:,1]))
