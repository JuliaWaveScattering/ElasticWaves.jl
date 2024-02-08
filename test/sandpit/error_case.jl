using ElasticWaves
using LinearAlgebra

#=
Create this code to test Tikanov with delta=0 in the test tha fails.
It fails at the frequencies 10 and 40, for basis order higher than 1 in the first case and higher than 2 in the second. 
Here I observed that the proporcionality of the fields  ||ϕ|| ∝||σ||/(ρ ω^2) breaks when we go to low frequency and high basis order. 
=#

#properties of the bearing
steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0)
ρ=steel.ρ


cs=Real(steel.cs)

#choosing basis_order

basis_order = 5

#minimum frequency for the diference in the radius

ω_min=max(Real((2*(basis_order)*(basis_order+1))^(1/2)*cs)/10, 2pi*cs/(bearing.inner_radius-bearing.outer_radius) )

#test ωs

ωs = [10.0,1*40,ω_min,1e4,4e4,7e4,6e7]

cs^2



# this non-dimensional number determines what basis_order is neeeded
kpas = bearing.outer_radius .* ωs ./ steel.cp
ksas = bearing.inner_radius .* ωs ./ steel.cs

#replicate of the test

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


#select the frequency from the vector of frequencies.

i=1

m=1

Ms = [
    boundarycondition_system(ωs[i], bearing, sims[i].boundarydata1.boundarytype, sims[i].boundarydata2.boundarytype, m)  
    for m in -basis_order:basis_order]

A = Ms[m]
b =[ [
    sims[i].boundarydata1.fourier_modes[m+basis_order+1,:];
    sims[i].boundarydata2.fourier_modes[m+basis_order+1,:] ] for m in -basis_order:basis_order
]

x = A \ b[m]

norm(x)

norm(b[m])/(ωs[i]^2*ρ)

n=basis_order

abs((norm(x)-norm(b[m])/(ωs[i]^2*ρ))/norm(x))


abs((norm(x)-norm(b[m])/(cs^2*n*(n+1)*ρ))/norm(x))


α=1/(ωs[i]^2*ρ)


δt=0.0

norm(A*x-b[m])^2


#δ=(0.2/(ksas[i])^2)


xT=[A; sqrt(δt) * I] \ [b[1]; zeros(size(A)[2])]

norm(x-xT)/norm(x)

 κ=cond(A)

 



wave = ElasticWave(sims[i])

ϕ=wave.potentials[1].coefficients

ψ=wave.potentials[2].coefficients

norm(ϕ[:,m])

norm(ψ[:,m])

norm(vcat(ϕ[:,m],ψ[:,m]))











