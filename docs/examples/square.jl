using ElasticWaves
using MultipleScattering
using Statistics
using Plots
using FFTW
using LinearAlgebra


basis_order =10;

numberofsensors = 12

basis_length = 2*basis_order + 1

ωs=[10.0,1e4,1e6]
#Friction coefficient

μ=0.5
#μ=1
θs = LinRange(0, 2pi, 10*(basis_length) + 1)[1:end-1]
#θs = LinRange(0, 2pi, 500)[1:end-1] 
#θ2s = LinRange(0, 2pi, 4*basis_length + 1)[1:end-1]
θ2s = LinRange(0, 2pi, 2000)[1:end-1] 
θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]

θ_line= pi

l=2

fourier_modes_p =[l *exp(-im*n*θ_line)* sinc(n*l/2pi)/2pi for n in -basis_order:basis_order ]

fourier_modes_s = μ.*[l *exp(-im*n*θ_line)* sinc(n*l/2pi)/2pi for n in -basis_order:basis_order ]

bc1_forward = TractionBoundary(inner=true)
bc2_forward = TractionBoundary(outer=true)

bc1_inverse = DisplacementBoundary(outer=true)
bc2_inverse = TractionBoundary(outer=true)

fouter= 0*fourier_modes_p

bd1_forward =  BoundaryData(bc1_forward, θs=θs, fourier_modes=hcat(fourier_modes_p,fourier_modes_s)) 
bd2_forward=BoundaryData(bc2_forward,θs=θs, fourier_modes=hcat(fouter,fouter))

i=2

sim = BearingSimulation(ωs[i], bearing, bd1_forward, bd2_forward; basis_order = basis_order)

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


scale = steel.ρ * ωs[i]^2

wave.potentials[1].coefficients

potential = HelmholtzPotential{2}(wave.potentials[1].wavespeed, wave.potentials[1].wavenumber, scale .* wave.potentials[1].coefficients)

res = field(potential, bearing; res = 120)

# plot the radial traction
plot(res,ωs[i]; seriestype=:heatmap, field_apply = f -> real(f[1]))
plot!(Circle(bearing.inner_radius))
plot!(Circle(bearing.outer_radius))

#INVERSE WITHOUT PRIOR

displacement_outer = [displacement(wave,x) for x in x_outer];
displacement_outer = hcat(displacement_outer...) |> transpose |> collect;


traction_outer= [traction(wave,x) for x in x_outer]
traction_outer = hcat(traction_outer...) |> transpose |> collect;


#Generate displacement for the inverse problem

bd1_inverse = BoundaryData(
    bc1_inverse;
    θs = θs_inv,
    fields = displacement_outer
)

bd1_inverse_modes=fields_to_fouriermodes(bd1_inverse, basis_order)
bd1_inverse_fields=fouriermodes_to_fields(bd1_inverse_modes)
norm(bd1_inverse.fields-bd1_inverse_fields.fields)

# traction_outer_forward= traction_outer_inverse

bd2_inverse= BoundaryData(
    bc2_inverse;
    θs = θs_inv,
    fields = traction_outer
)

bd2_inverse_modes=fields_to_fouriermodes(bd2_inverse, basis_order)
bd2_inverse_fields=fouriermodes_to_fields(bd2_inverse_modes)
norm(bd2_inverse.fields-bd2_inverse_fields.fields)

#bd2_inverse=fields_to_fouriermodes(bd2_inverse, basis_order)

inverse_sim = BearingSimulation(ωs[i], bearing, bd1_inverse, bd2_inverse, basis_order=basis_order)    
 res = field(wave, bearing, TractionType(); res = 70)

inv_wave=ElasticWave(inverse_sim)

#Calculate and compute for the inverse problem without prior.

traction_for_ω= [traction(wave,x) for x in x2_inner]
traction_for_ω = hcat(traction_for_ω...) |> transpose |> collect

traction_inv_ω = [traction(inv_wave,x) for x in x2_inner];
traction_inv_ω = hcat(traction_inv_ω...) |> transpose |> collect


plot(θ2s,real.(traction_for_ω[:,1]), linestyle = :dash, linewidth = 2)
plot!(θ2s,real.(traction_inv_ω[:,1]), linestyle = :dash, linewidth = 2)


plot(θ2s,real.(traction_for_ω[:,2]), linestyle = :dash, linewidth = 2)
plot!(θ2s,real.(traction_inv_ω[:,2]), linestyle = :dash, linewidth = 2)



scale = steel.ρ * ωs[i]^2

inv_wave.potentials[1].coefficients

potential = HelmholtzPotential{2}(inv_wave.potentials[1].wavespeed, 
inv_wave.potentials[1].wavenumber, scale .* inv_wave.potentials[1].coefficients)

res = field(potential, bearing; res = 120)

plot(res,ωs[i]; seriestype=:heatmap, field_apply = f -> real(f[1]))
plot!(Circle(bearing.inner_radius))
plot!(Circle(bearing.outer_radius))

#PRIOR METHOD



Fp=[l *exp(-im*n*θ_line)* sinc(n*l/2pi)/2pi for n in -basis_order:basis_order ] 
Fs=[l *exp(-im*n*θ_line)* sinc(n*l/2pi)/2pi for n in -basis_order:basis_order ]




f0=0.0*Fp




bd1= BoundaryData(bc1_forward, θs=θs, fourier_modes=hcat(Fp,f0)) 
bd2= BoundaryData(bc1_forward, θs=θs, fourier_modes=hcat(f0,Fs)) 



boundarybasis=  BoundaryBasis( [ bd1 , bd2] )   

#Generate boundary data

bd1_inverse = BoundaryData(
    bc1_inverse;
    θs = θs_inv,
    fields = displacement_outer
)


#bd1_inverse=fields_to_fouriermodes(bd1_inverse, basis_order)


bd2_inverse= BoundaryData(
    bc2_inverse;
    θs = θs_inv,
    fields = traction_outer
)

#bd2_inverse=fields_to_fouriermodes(bd2_inverse, basis_order)

inverse_sim = BearingSimulation(ωs[i], bearing, bd1_inverse, bd2_inverse, boundarybasis1=boundarybasis,basis_order=basis_order)    
# res = field(wave, bearing, TractionType(); res = 70)

inv_wave=ElasticWave(inverse_sim)

#solve with prior.

#calculate_predicted_tractions

traction_inv_ω = [traction(inv_wave,x) for x in x2_inner];
traction_inv_ω = hcat(traction_inv_ω...) |> transpose |> collect


plot(θ2s,real.(traction_for_ω[:,1]), linestyle = :dash, linewidth = 2)
plot!(θ2s,real.(traction_inv_ω[:,1]), linestyle = :dash, linewidth = 2)

plot(θ2s,real.(traction_for_ω[:,2]), linestyle = :dash, linewidth = 2)
plot!(θ2s,real.(traction_inv_ω[:,2]), linestyle = :dash, linewidth = 2)

maximum(real.(traction_inv_ω[:,2]))/maximum(real.(traction_inv_ω[:,1]))

#plot the fields

scale = steel.ρ * ωs[i]^2


inv_wave.potentials[1].coefficients

potential = HelmholtzPotential{2}(inv_wave.potentials[1].wavespeed, 
inv_wave.potentials[1].wavenumber, scale .* inv_wave.potentials[1].coefficients)

res = field(potential, bearing; res = 120)

plot(res,ωs[i]; seriestype=:heatmap, field_apply = f -> real(f[1]))
plot!(Circle(bearing.inner_radius))
plot!(Circle(bearing.outer_radius))
