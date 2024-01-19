using ElasticWaves
using MultipleScattering
using Statistics
using Plots
using FFTW
using LinearAlgebra

ω=1e6
basis_order = 10;
numberofsensors = 13
basis_length = 2*basis_order + 1
μ=0.081
#μ=1
#θs = LinRange(0, 2pi, basis_length + 1)[1:end-1]
θs = LinRange(0, 2pi, 500)[1:end-1] 
#θ2s = LinRange(0, 2pi, 4*basis_length + 1)[1:end-1]
θ2s = LinRange(0, 2pi, 2000)[1:end-1] 
θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]

#number_of_ωs=500

ωs=LinRange(10,1e5, 100)

steel = Elastic(2; ρ = 7.80, cp = 5.00, cs = 3.50)
bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0, number_of_rollers=1.0)
Z=bearing.number_of_rollers
Δt = real((bearing.outer_radius - bearing.inner_radius)/steel.cp)
d=2pi/Z
Ω=10.0
tmax =1*d/Ω
#tmax=Δt
#tmax=0.05
number_of_ts= 1*(tmax/Δt)
number_of_ts= Int(length(taux))


#t=LinRange(0,150, number_of_ts)|>transpose
#taux=LinRange(0,tmax, Int(number_of_ts))|>transpose |>collect
taux=ω_to_t(ωs) |>transpose |>collect
t=repeat(taux, outer= length(θs))
#Ω=2pi/tmax


n_order=10
# the pressure and shear fields 
#fp1 = 1*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im
#fs1 = 1*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im


#FORWARD PROBLEM
fp= [(Z/(2pi*Ω)).* sum(exp.(im.*n.*Z.*(Ω.*t[:,i].-θs)) for n in -n_order:n_order ) for i in 1:number_of_ts ]
fs=μ.* fp
#fp=0*fp
anim= @animate for i in 1:number_of_ts
    plot(θs,real.(fp[i]),legend=false, ylims=(-0.1,0.4),xlims=(0,2pi))
    plot!(θs,real.(fs[i]),legend=false, ylims=(-0.1,0.4),xlims=(0,2pi))
    

end

gif(anim,"docs/images/delta_forcing.gif", fps = 10)

fp=hcat(fp...)|> transpose 
fs=hcat(fs...)|> transpose

taux |> transpose |>collect
taux=LinRange(0,tmax, Int(number_of_ts))
Fp=time_to_frequency(fp,taux)
Fs=μ.*Fp

#Fp=0*Fp
#ωs0 = t_to_ω(taux)
#ωs = LinRange(0.02,maximum(ωs0),length(ωs0))
bc1_forward = TractionBoundary(inner=true)
bc2_forward = TractionBoundary(outer=true)

bc1_inverse = DisplacementBoundary(outer=true)
bc2_inverse = TractionBoundary(outer=true)

fouter= 0*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

bd1_forward=[ BoundaryData(bc1_forward, θs=θs, fields=hcat(Fp[i,:],Fs[i,:])) for i in 1:length(ωs) ]

bd2_forward= BoundaryData(bc2_forward,θs=θs, fields=hcat(fouter,fouter))

#i=175
k=1
indexes= [findall(x->x==k*Ω, ωs) for k in 1:5]

i=1

ωs[i]

plot(θs, abs.(Fp[i,:]))
plot(θs, abs.(Fs[i,:]))


sim = BearingSimulation(ωs[i], bearing, bd1_forward[i], bd2_forward;
basis_order = basis_order)
wave = ElasticWave(sim)

# res = field(wave, bearing, TractionType(); res = 70)

# scale the potential to match the units of stress
scale = steel.ρ * ωs[i]^2

wave.potentials[1].coefficients

potential = HelmholtzPotential{2}(wave.potentials[1].wavespeed, wave.potentials[1].wavenumber, scale .* wave.potentials[1].coefficients)

res = field(potential, bearing; res = 120)

# plot the radial traction
plot(res,ωs[i]; seriestype=:heatmap, field_apply = f -> real(f[1]))
plot!(Circle(bearing.inner_radius))
plot!(Circle(bearing.outer_radius))


x_outer=[radial_to_cartesian_coordinates([bearing.outer_radius,θ]) for θ in θs_inv ]

x2_inner = [
    radial_to_cartesian_coordinates([bearing.inner_radius, θ])
for θ in θ2s]


x2_outer = [
    radial_to_cartesian_coordinates([bearing.outer_radius, θ])
for θ in θ2s]


displacement_outer = [displacement(wave,x) for x in x_outer];
displacement_outer = hcat(displacement_outer...) |> transpose |> collect;



bd1_inverse = BoundaryData(
bc1_inverse;
θs = θs_inv,
fields = displacement_outer
)

#bd1_inverse=fields_to_fouriermodes(bd1_inverse, basis_order)

bd2_inverse= bd2_forward

inverse_sim = BearingSimulation(ωs[i], bearing, bd1_inverse, bd2_inverse)    
# res = field(wave, bearing, TractionType(); res = 70)

inv_wave=ElasticWave(inverse_sim)

traction_inv_ω = [traction(inv_wave,x) for x in x2_inner];
traction_inv_ω = hcat(traction_inv_ω...) |> transpose |> collect

plot(θs, abs.(Fp[i,:]))
plot!(θ2s,real.(traction_inv_ω[:,1]))



plot(θs, abs.(Fs[i,:]))
plot!(θ2s,real.(traction_inv_ω[:,2]))


scale = steel.ρ * ωs[i]^2

inv_wave.potentials[1].coefficients

potential = HelmholtzPotential{2}(inv_wave.potentials[1].wavespeed, 
inv_wave.potentials[1].wavenumber, scale .* inv_wave.potentials[1].coefficients)

res = field(potential, bearing; res = 120)

plot(res,ωs[i]; seriestype=:heatmap, field_apply = f -> real(f[1]))
plot!(Circle(bearing.inner_radius))
plot!(Circle(bearing.outer_radius))


#PRIOR METHOD

fp1= [(Z/(2pi*Ω)).* sum(exp.(im.*n.*Z.*(Ω.*t[:,i].-θs)) for n in -n_order:n_order ) for i in 1:number_of_ts ]
fs1= fp1

fp1=hcat(fp1...)|> transpose 
fs1=hcat(fs1...)|> transpose

Fp1=time_to_frequency(fp1,taux)
Fs1=time_to_frequency(fs1,taux)

f0=0 .*Fp1

#Creating the basis for each frequency

bd1=[ BoundaryData(bc1_forward, θs=θs, fields=hcat(Fp1[i,:],f0[i,:])) for i in 1:length(ωs) ]
bd2=[ BoundaryData(bc1_forward, θs=θs, fields=hcat(f0[i,:],Fs1[i,:])) for i in 1:length(ωs) ]
#plot(θs,abs.(Fp1[50,:]))

boundarybasis= [ BoundaryBasis( [ bd1[i] , bd2[i]] ) for i in 1:length(ωs)  ]

bd1_inverse = BoundaryData(
    bc1_inverse;
    θs = θs_inv,
    fields = displacement_outer
)

#bd1_inverse=fields_to_fouriermodes(bd1_inverse, basis_order)

bd2_inverse= bd2_forward

inverse_sim = BearingSimulation(ωs[i], bearing, bd1_inverse, bd2_inverse, boundarybasis1=boundarybasis[i])    
# res = field(wave, bearing, TractionType(); res = 70)

inv_wave=ElasticWave(inverse_sim)

traction_inv_ω = [traction(inv_wave,x) for x in x2_inner];
traction_inv_ω = hcat(traction_inv_ω...) |> transpose |> collect

plot(θs, abs.(Fp[i,:]))
plot!(θ2s,real.(traction_inv_ω[:,1]))



plot(θs, abs.(Fs[i,:]))
plot!(θ2s,real.(traction_inv_ω[:,2]))


scale = steel.ρ * ωs[i]^2

inv_wave.potentials[1].coefficients

potential = HelmholtzPotential{2}(inv_wave.potentials[1].wavespeed, 
inv_wave.potentials[1].wavenumber, scale .* inv_wave.potentials[1].coefficients)

res = field(potential, bearing; res = 120)

plot(res,ωs[i]; seriestype=:heatmap, field_apply = f -> real(f[1]))
plot!(Circle(bearing.inner_radius))
plot!(Circle(bearing.outer_radius))
