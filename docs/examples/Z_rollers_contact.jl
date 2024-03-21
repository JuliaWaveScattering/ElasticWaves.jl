
using ElasticWaves
using MultipleScattering
using Statistics
using Plots
using FFTW
using LinearAlgebra

#=
(For 1 Roller)
The inverse problem solves well (almost exact) for number_of_sensors= basis_length for both cases (with prior and without). 
If we put less than that the tractions observed in the inverse problem without prior are terrible. 
The prior method can give good results with less than that. with number_of_sensors= basis_order+2 it gives good results
 (For 10 Rollers)
 with basis_order=10 and number_of_sensors=21=basis_length we can see the deltas but not as in the input.
=#

#Definining parameters of the problem
#ω=1e6
basis_order = 10;
numberofsensors = 3
basis_length = 2*basis_order + 1

#Friction coefficient

μ=0.081
μ=0.5
#θs = LinRange(0, 2pi, basis_length + 1)[1:end-1]
θs = LinRange(0, 2pi, 500)[1:end-1] 
#θ2s = LinRange(0, 2pi, 4*basis_length + 1)[1:end-1]
θ2s = LinRange(0, 2pi, 2000)[1:end-1] 
θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]

#Properties of the bearing
steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0, number_of_rollers=1.0)
Z=bearing.number_of_rollers
Δt = real((bearing.outer_radius - bearing.inner_radius)/steel.cp)
d=2pi/Z

#Angular velocity
Ω=10

#time that we are observing the waves (one period of the bearings rolling)

tmax =1*d/Ω
tmax=Δt
#tmax=0.05
number_of_ts= 1*(tmax/Δt)
number_of_ts= 124


#t=LinRange(0,150, number_of_ts)|>transpose
taux=LinRange(0,tmax, Int(number_of_ts))|>transpose |>collect
t=repeat(taux, outer= length(θs))
#Ω=2pi/tmax


n_order=basis_order
# the pressure and shear fields 
#fp1 = 1*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im
#fs1 = 1*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im


#FORWARD PROBLEM
#Construct the forcing fp(θ,t)= Z/(2πΩ) ∑_n  exp(iZn(Ωt-θ) )

fp= [(Z/(2pi*Ω)).* sum(exp.(im.*n.*Z.*(Ω.*t[:,i].-θs)) for n in -n_order:n_order ) for i in 1:number_of_ts ]  
fs=μ.* fp
#fp=0*fp
#gif of the forcing
anim= @animate for i in 1:number_of_ts
    plot(θs,real.(fp[i]),legend=false, ylims=(-0.1,0.4),xlims=(0,2pi))
    plot!(θs,real.(fs[i]),legend=false, ylims=(-0.1,0.4),xlims=(0,2pi))
    

end

gif(anim,"docs/images/delta_forcing.gif", fps = 10)

# fourier transform of the force
fp=hcat(fp...)|> transpose 
fs=hcat(fs...)|> transpose

taux |> transpose |>collect
taux=LinRange(0,tmax, Int(number_of_ts))
Fp=time_to_frequency(fp,taux)
Fs=μ.*Fp

#Fp=0*Fp

#Get frequencies
ωs0 = t_to_ω(taux)
#Avoiding ω=0
ωs = LinRange(0.02,maximum(ωs0),length(ωs0))

#Setting Boundary conditions

bc1_forward = TractionBoundary(inner=true)
bc2_forward = TractionBoundary(outer=true)

bc1_inverse = DisplacementBoundary(outer=true)
bc2_inverse = TractionBoundary(outer=true)

#FORWARD PROBLEM

fouter= 0*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

bd1_forward=[ BoundaryData(bc1_forward, θs=θs, fields=hcat(Fp[i,:],Fs[i,:])) for i in 1:length(ωs0) ]

bd2_forward= BoundaryData(bc2_forward,θs=θs, fields=hcat(fouter,fouter))

#solving for each frequency

results = map(eachindex(ωs)) do i
    sim = BearingSimulation(ωs[i], bearing, bd1_forward[i], bd2_forward;
        basis_order = basis_order)
    wave = ElasticWave(sim)

    # res = field(wave, bearing, TractionType(); res = 70)

    # scale the potential to match the units of stress
    scale = steel.ρ * ωs[i]^2

    wave.potentials[1].coefficients

    potential = HelmholtzPotential{2}(wave.potentials[1].wavespeed, wave.potentials[1].wavenumber, scale .* wave.potentials[1].coefficients, wave.potentials[1].modes)

    res = field(potential, bearing; res = 120)
    # plot(res, ωs[i]; seriestype=:coω_to_t(ωs)ntour)
end


#calculate the fields

fields = [[f[1] for f in field(r)] for r in results];
all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs);




Δt = (bearing.outer_radius - bearing.inner_radius) / steel.cp |> real

# need at as many samples as possible within this time frame
dt = Δt / 4

# from dt we can calculate what's the maximum frequency we need
maxω = 2pi / dt
impulse = GaussianImpulse(maxω)
ts = ω_to_t(ωs)
time_result = frequency_to_time(all_results; t_vec = ts,impulse=impulse)

maxc = 0.33 .* maximum(norm.(field(time_result)))
minc = - maxc

minimum(norm.(field(time_result)))

#fp=10*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im
#fs=5 *exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

r = bearing.inner_radius / 4

    anim = @animate for t in ts
        plot(time_result, t,
          seriestype=:heatmap,
          clim = (1*minc, 1*maxc),
          leg = false,
        )
        plot!(frame = :none, title="", xguide ="",yguide ="")
        plot!(Circle(bearing.inner_radius))
        plot!(Circle(bearing.outer_radius))
        #plot!(Circle([r-bearing.inner_radius, 0.0], r))
        plot!(Circle(bearing.inner_radius -2r))
    end
    
    gif(anim,"docs/images/pressure-Zrollers.gif", fps = 10)

#gif of the fields

#INVERSE PROBLEM WITHOUT PRIOR


x_outer=[radial_to_cartesian_coordinates([bearing.outer_radius,θ]) for θ in θs_inv ]

#x2_inner are the cooordinates that we will calculate the forces that the fields of our solutions generates

x2_inner = [
    radial_to_cartesian_coordinates([bearing.inner_radius, θ])
for θ in θ2s]


x2_outer = [
    radial_to_cartesian_coordinates([bearing.outer_radius, θ])
for θ in θ2s]


#vector to store tractions for each frequency
traction_inv=[]

#traction_inv_out = []

#solve the inverse problem for each frequency

results = map(eachindex(ωs)) do i
    sim = BearingSimulation(ωs[i], bearing, bd1_forward[i], bd2_forward;
        basis_order = basis_order)
    wave = ElasticWave(sim)

    displacement_outer = [displacement(wave,x) for x in x_outer];
    displacement_outer = hcat(displacement_outer...) |> transpose |> collect;

    

    bd1_inverse = BoundaryData(
        bc1_inverse;
        θs = θs_inv,
        fields = displacement_outer
    )

    #bd1_inverse=fields_to_fouriermodes(bd1_inverse, basis_order)

    bd2_inverse= bd2_forward

    inverse_sim = BearingSimulation(ωs[i], bearing, bd1_inverse, bd2_inverse, basis_order=basis_order)    
    # res = field(wave, bearing, TractionType(); res = 70)

    inv_wave=ElasticWave(inverse_sim)
    
    #calculate traction for each frequency
    traction_inv_ω = [traction(inv_wave,x) for x in x2_inner];
    traction_inv_ω = hcat(traction_inv_ω...) |> transpose |> collect

    #save each traction for each frequency
    push!(traction_inv,traction_inv_ω )

    #traction_inv_out_aux = [traction(inv_wave,x) for x in x2_outer];
    #traction_inv_out_aux = hcat(traction_inv_out_aux...) |> transpose |> collect


    # scale the potential to match the units of stress
    scale = steel.ρ * ωs[i]^2

    inv_wave.potentials[1].coefficients

    potential = HelmholtzPotential{2}(inv_wave.potentials[1].wavespeed, 
    inv_wave.potentials[1].wavenumber, scale .* inv_wave.potentials[1].coefficients, inv_wave.potentials[1].modes)

    res = field(potential, bearing; res = 120)
    # plot(res, ωs[i]; seriestype=:coω_to_t(ωs)ntour)
end


#seperate pressure and shear forces for each frequency

traction_inv_p=[]

traction_inv_s=[]

for i in 1:length(ωs)

    traction_inv_p_ω=traction_inv[i][:,1]

    push!(traction_inv_p, traction_inv_p_ω)

    traction_inv_s_ω=traction_inv[i][:,2]

    push!(traction_inv_s, traction_inv_s_ω)

end

traction_inv_p= hcat(traction_inv_p...) |> transpose |>collect

traction_inv_s= hcat(traction_inv_s...) |> transpose |>collect

#fourier transform back to time
fp_inv= frequency_to_time(traction_inv_p, ωs)

fs_inv= frequency_to_time(traction_inv_s, ωs)

maxf=maximum(fp_inv)
minf=minimum(fs_inv)

# gif of the forces
anim= @animate for i in 1:length(ts)
    plot(θ2s,real.(fp_inv[i,:]),legend=false, ylims=(-0.1,0.4),xlims=(0,2pi))
    plot!(θ2s,real.(fs_inv[i,:]),legend=false, ylims=(-0.1,0.4),xlims=(0,2pi))
    

end


gif(anim,"docs/images/delta_forcing_inv.gif", fps = 10)


# doing gif of the fields

fields = [[f[1] for f in field(r)] for r in results];
all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs);




Δt = (bearing.outer_radius - bearing.inner_radius) / steel.cp |> real

# need at as many samples as possible within this time frame
dt = Δt / 4

# from dt we can calculate what's the maximum frequency we need
maxω = 2pi / dt
impulse = GaussianImpulse(maxω)
ts = ω_to_t(ωs)
time_result = frequency_to_time(all_results; t_vec = ts, impulse=impulse)

maxc = 0.33 .* maximum(norm.(field(time_result)))
minc = - maxc

minimum(norm.(field(time_result)))

#fp=10*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im
#fs=5 *exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

r = bearing.inner_radius / 4

    anim = @animate for t in ts
        plot(time_result, t,
          seriestype=:heatmap,
          clim = (1*minc, 1*maxc),
          leg = false,
        )
        plot!(frame = :none, title="", xguide ="",yguide ="")
        plot!(Circle(bearing.inner_radius))
        plot!(Circle(bearing.outer_radius))
        #plot!(Circle([r-bearing.inner_radius, 0.0], r))
        plot!(Circle(bearing.inner_radius -2r))
    end
    
    gif(anim,"docs/images/pressure-Zrollers-inverse.gif", fps = 10)



#INVERSE PROBLEM WITH PRIOR

#contructing prior, I am saying to my code that the force is a sequence of deltas walking but not giving the amplitudes.

fp1= [(Z/(2pi*Ω)).* sum(exp.(im.*n.*Z.*(Ω.*t[:,i].-θs)) for n in -n_order:n_order ) for i in 1:number_of_ts ]
fs1= fp1

fp1=hcat(fp1...)|> transpose 
fs1=hcat(fs1...)|> transpose

#Writing my prior in frequency space
Fp1=time_to_frequency(fp1,taux)
Fs1=time_to_frequency(fs1,taux)

#Creating a vector of zeros
f0=0 .*Fp1

#Creating the basis for each frequency, first crate the boundary data of the prior.

bd1=[ BoundaryData(bc1_forward, θs=θs, fields=hcat(Fp1[i,:],f0[i,:])) for i in 1:length(ωs) ]
bd2=[ BoundaryData(bc1_forward, θs=θs, fields=hcat(f0[i,:],Fs1[i,:])) for i in 1:length(ωs) ]

#plot(θs,abs.(Fp1[50,:]))

#Create the boundary_basis

boundarybasis= [ BoundaryBasis( [ bd1[i] , bd2[i]] ) for i in 1:length(ωs)  ]

#from here is the same path of the previous case but passing the boundary_basis

x_outer=[radial_to_cartesian_coordinates([bearing.outer_radius,θ]) for θ in θs_inv ]

x2_inner = [
    radial_to_cartesian_coordinates([bearing.inner_radius, θ])
for θ in θ2s]


x2_outer = [
    radial_to_cartesian_coordinates([bearing.outer_radius, θ])
for θ in θ2s]


traction_inv=[]

#traction_inv_out = []



results = map(eachindex(ωs)) do i
    sim = BearingSimulation(ωs[i], bearing, bd1_forward[i], bd2_forward;
        basis_order = basis_order)
    wave = ElasticWave(sim)

    displacement_outer = [displacement(wave,x) for x in x_outer];
    displacement_outer = hcat(displacement_outer...) |> transpose |> collect;

    

    bd1_inverse = BoundaryData(
        bc1_inverse;
        θs = θs_inv,
        fields = displacement_outer
    )

    #bd1_inverse=fields_to_fouriermodes(bd1_inverse, basis_order)

    bd2_inverse= bd2_forward

    inverse_sim = BearingSimulation(ωs[i], bearing, bd1_inverse, bd2_inverse, boundarybasis1=boundarybasis[i];
     basis_order=basis_order)    
    # res = field(wave, bearing, TractionType(); res = 70)

    inv_wave=ElasticWave(inverse_sim)

    traction_inv_ω = [traction(inv_wave,x) for x in x2_inner];
    traction_inv_ω = hcat(traction_inv_ω...) |> transpose |> collect

    push!(traction_inv,traction_inv_ω )

    #traction_inv_out_aux = [traction(inv_wave,x) for x in x2_outer];
    #traction_inv_out_aux = hcat(traction_inv_out_aux...) |> transpose |> collect


    # scale the potential to match the units of stress
    scale = steel.ρ * ωs[i]^2

    inv_wave.potentials[1].coefficients

    potential = HelmholtzPotential{2}(inv_wave.potentials[1].wavespeed, 
    inv_wave.potentials[1].wavenumber, scale .* inv_wave.potentials[1].coefficients, inv_wave.potentials[1].modes)

    res = field(potential, bearing; res = 120)
    # plot(res, ωs[i]; seriestype=:coω_to_t(ωs)ntour)
end


traction_inv_p=[]

traction_inv_s=[]

for i in 1:length(ωs)

    traction_inv_p_ω=traction_inv[i][:,1]

    push!(traction_inv_p, traction_inv_p_ω)

    traction_inv_s_ω=traction_inv[i][:,2]

    push!(traction_inv_s, traction_inv_s_ω)

end

traction_inv_p= hcat(traction_inv_p...) |> transpose |>collect

traction_inv_s= hcat(traction_inv_s...) |> transpose |>collect


fp_inv= frequency_to_time(traction_inv_p, ωs)

fs_inv= frequency_to_time(traction_inv_s, ωs)

maxf=maximum(fp_inv)
minf=minimum(fs_inv)

anim= @animate for i in 1:length(ts)
    plot(θ2s,real.(fp_inv[i,:]),legend=false, ylims=(-0.1,0.4),xlims=(0,2pi))
    plot!(θ2s,real.(fs_inv[i,:]),legend=false, ylims=(-0.1,0.4),xlims=(0,2pi))
    

end

plot(θ2s,real.(fp_inv[13,:]),legend=false,xlims=(0,2pi))


gif(anim,"docs/images/delta_forcing_inv_prior.gif", fps = 10)




fields = [[f[1] for f in field(r)] for r in results];
all_results = FrequencySimulationResult(hcat(fields...), results[1].x, ωs);




Δt = (bearing.outer_radius - bearing.inner_radius) / steel.cp |> real

# need at as many samples as possible within this time frame
dt = Δt / 4

# from dt we can calculate what's the maximum frequency we need
maxω = 2pi / dt
impulse = GaussianImpulse(maxω)
ts = ω_to_t(ωs)
time_result = frequency_to_time(all_results; t_vec = ts, impulse=impulse)

maxc = 0.33 .* maximum(norm.(field(time_result)))
minc = - maxc

minimum(norm.(field(time_result)))

#fp=10*exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im
#fs=5 *exp.(-20.0 .* (θs .- pi).^2) + θs .* 0im

r = bearing.inner_radius / 4

    anim = @animate for t in ts
        plot(time_result, t,
          seriestype=:heatmap,
          clim = (1*minc, 1*maxc),
          leg = false,
        )
        plot!(frame = :none, title="", xguide ="",yguide ="")
        plot!(Circle(bearing.inner_radius))
        plot!(Circle(bearing.outer_radius))
        #plot!(Circle([r-bearing.inner_radius, 0.0], r))
        plot!(Circle(bearing.inner_radius -2r))
    end
    
    gif(anim,"docs/images/pressure-Zrollers-inverse-prior.gif", fps = 10)

