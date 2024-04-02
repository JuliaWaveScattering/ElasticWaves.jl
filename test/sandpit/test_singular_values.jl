using ElasticWaves
# using Plots
using LinearAlgebra
using BlockDiagonals

steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)

steel = Elastic(2; ρ = 1.0, cp = 3.0, cs = 2.0)
bearing = RollerBearing(medium = steel, inner_radius = 1.5, outer_radius = 2.0)

bc_inv1 = DisplacementBoundary(outer=true)
bc_inv2 = TractionBoundary(outer=true)
bc_for1 = TractionBoundary(inner=true)
bc_for2= TractionBoundary(outer=true)

ω = 0.5
n = 20;

cp = steel.cp; cs = steel.cs
kP = ω / cp; kS = ω / cs

M = boundarycondition_system(ω, bearing, bc_for1, bc_for2, n)

cond(M)

r1 = bearing.inner_radius;
r2 = bearing.outer_radius;

S = diagm([
    besselj(n, kP * r1)/2 + besselj(n, kP * r2)/2, 
    hankelh1(n, kP * r1)/2 + hankelh1(n, kP * r2)/2,
    besselj(n, kS * r1)/2 + besselj(n, kS * r2)/2, 
    hankelh1(n, kS * r1)/2 + hankelh1(n, kS * r2)/2,
])

cond(M * inv(S))

ωs = LinRange(1.0e1,1.0e6,1000)
basis_order = 5
basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

kpas = ωs .*(bearing.outer_radius-bearing.inner_radius) ./ steel.cp |>real

M_inv = [boundarycondition_system(ω, bearing, bc_inv1, bc_inv2, i) for ω in ωs, i in -basis_order:basis_order]

M_inv[1,:]

Minv=[BlockDiagonal(M_inv[i,:]) for i in 1:length(ωs)]


#svdMs=svd.(M_inv)

svdvalues=svdvals.(Minv)

cond_numbers= cond.(Minv)

#svdvals_min=zeros(length(ωs), basis_length)

plot(ωs,minimum.(svdvalues), yaxis=:log, xaxis=:log )
xlabel!("ω")
ylabel!("minimum singular value")


plot(ωs,cond_numbers, yaxis=:log, xaxis=:log )
xlabel!("ω")
ylabel!("conditioning number")


#=
for i in  1:length(ωs)
    for j in  1:basis_length
        svdvals_min[i,j]=minimum(svdMs[i,j].S)
    end
end

plot(ωs,svdvals_min[:,1],xaxis=:log,yaxis=:log)
xlabel!("ω")
ylabel!("minimum singular value per mode")

plot!(ωs,svdvals_min[:,2])

plot!(ωs,svdvals_min[:,3])
plot!(ωs,svdvals_min[:,4])
plot!(ωs,svdvals_min[:,5])
plot!(ωs,svdvals_min[:,6])
=#

forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im

bd1= BoundaryData(TractionBoundary(inner = true); coefficients =  hcat(forcing_modes[:,1],0*forcing_modes[:,1]))
bd2= BoundaryData(TractionBoundary(inner = true); coefficients =  hcat(0*forcing_modes[:,2],forcing_modes[:,2]))

boundarybasis1=BoundaryBasis([bd1,bd2])

bd3= BoundaryData(TractionBoundary(outer  = true); coefficients =  hcat(forcing_modes[:,3],0*forcing_modes[:,3]))
bd4= BoundaryData(TractionBoundary(outer  = true); coefficients =  hcat(0*forcing_modes[:,4],forcing_modes[:,4]))

boundarybasis2=BoundaryBasis([bd3,bd4])
    

    # Create the Prior matrix P 
N1 = length(boundarybasis1.basis)
F1s = [b.coefficients for b in boundarybasis1.basis]

F1 = vcat((transpose.(F1s))...)
P1s = [reshape(F1[:,m],2,N1) for m = 1:size(F1,2)]
P1 = vcat([ [P; 0.0*P]  for P in P1s]...)

N2 = length(boundarybasis2.basis)
F2s = [b.coefficients for b in boundarybasis2.basis]

F2 = vcat((transpose.(F2s))...)
P2s = [reshape(F2[:,m],2,N2) for m = 1:size(F2,2)]
P2 = vcat([ [0.0*P; P]  for P in P2s]...)

        # reshape just in case matrix is empty so hcat can work seemlesly
P1 = reshape(P1, 4 * (2 * basis_order + 1), :)
P2 = reshape(P2, 4 * (2 * basis_order + 1), :)
    
P = hcat(P1,P2)
    


M_for=[boundarycondition_system(ω, bearing, bc_for1, bc_for2, i) for ω in ωs, i in -basis_order:basis_order]

Mfor=[BlockDiagonal(M_for[i,:]) for i in 1:length(ωs)]

B=[Mfor[i]\P for i in 1:length(ωs)]

A=[Minv[i]*B[i] for i  in 1:length(ωs) ]

svdvalues_prior= svdvals.(A)

cond_numbers_prior=cond.(A)

plot(ωs,minimum.(svdvalues_prior), yaxis=:log, xaxis=:log )
xlabel!("ω")
ylabel!("minimum singular value")



plot(ωs,cond_numbers_prior, yaxis=:log, xaxis=:log )
xlabel!("ω")
ylabel!("conditioning number")
