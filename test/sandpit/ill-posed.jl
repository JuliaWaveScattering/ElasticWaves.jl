include("../../src/ElasticWaves.jl")

# choose a low frequency so it is ill-posed
    ω = 10.0

    steel = Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    bearing = RollerBearing(medium=steel, inner_radius=1.0, outer_radius = 2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpas = bearing.outer_radius * ω / steel.cp
    ksas = bearing.inner_radius * ω / steel.cs

    basis_order = 10
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)
 
## Stability check by adding Gaussian noise

    forcing_modes = rand(basis_length,4) + rand(basis_length,4) .* im
    bd1 = BoundaryData(TractionBoundary(inner=true); fourier_modes=forcing_modes[:, 1:2])
    bd2 = BoundaryData(TractionBoundary(outer=true); fourier_modes=forcing_modes[:, 3:4])

    sim = BearingSimulation(ω, bearing, bd1, bd2) 
    wave = ElasticWave(sim)
    
    # add 1% error to boundary conditions
    ε = 0.01 * maximum(abs.(forcing_modes));

    bd1.fourier_modes[:,:] = bd1.fourier_modes[:,:] + ε .* rand(basis_length, 2) + ε .* rand(basis_length, 2) .* im
    bd2.fourier_modes[:,:] = bd2.fourier_modes[:,:] + ε .* rand(basis_length, 2) + ε .* rand(basis_length, 2) .* im

    sim = BearingSimulation(ω, bearing, bd1, bd2) 
    waveδ = ElasticWave(sim)

    # I suppose 20% error is quite large. Any error in these coefficients will translate to an error in the field
    maximum(abs.(wave.potentials[1].coefficients - waveδ.potentials[1].coefficients)) / mean(abs.(wave.potentials[1].coefficients)) 
    maximum(abs.(wave.potentials[2].coefficients - waveδ.potentials[2].coefficients) ./ abs.(wave.potentials[2].coefficients))

    # predict the forcing modes  
# errors = map(eachindex(ωs)) do i

#     fs = map(-basis_order:basis_order) do m
#         coes = vcat(
#             waves[i].potentials[1].coefficients[:, m+basis_order+1],
#             waves[i].potentials[2].coefficients[:, m+basis_order+1]
#         )
#         boundarycondition_system(ωs[i], bearing, bd1.boundarytype, bd2.boundarytype, m) * coes
#     end
#     f = hcat(fs...) |> transpose
#     # maximum(abs.(f - hcat(bd1.fourier_modes,bd2.fourier_modes)))
#     maximum(abs.(f - forcing_modes))
# end

    f = hcat(
        field_modes(waveδ, bearing.inner_radius, bd1.boundarytype.fieldtype),
        field_modes(waveδ, bearing.outer_radius, bd1.boundarytype.fieldtype)
    )

    # gives great results when compared with the forcing modes used
    maximum(abs.(f - hcat(bd1.fourier_modes,bd2.fourier_modes)))

    # gives, as it should, poorer results when compared with the original forcing 
    maximum(abs.(f - forcing_modes))

    # Yet these result still seem too good. Let us look at the condition numbers for the systems we are solving
    σs = map(0:basis_order) do m
        A = boundarycondition_system(ω, bearing, bd1.boundarytype, bd2.boundarytype, m)
        cond(A)
    end

    # THe last one is larger than 1e50 ! Meaning that to solve A * x = b, an error of 1e-10 in b could result in an error of about  1e50 * 1e-10 = 1e40 when solving for x!

    # Let us see this in practice
    m = basis_order 
    A = boundarycondition_system(ω, bearing, bd1.boundarytype, bd2.boundarytype, m)

    x = rand(4)
    b = A * x

    # a direct solve gives a terrible solution, as it should.  
    x2 = A \ b
    norm(A * x2 - b)
 
    # a regulariser improves it, but nothing can fix this
    δ = 1e-8;
    Aδ = [A; I .* δ]
    bδ = [b; [0,0,0,0]]

    xδ = (Aδ \ bδ)
    
    norm(x - xδ)

    # So why did the method seem to work better than this? Let's us the same forcing from the bearing problem 

    b = hcat(bd1.fourier_modes, bd2.fourier_modes)[m + 1 + basis_order,:]
    x = [wave.potentials[1].coefficients[:, m+basis_order+1]; wave.potentials[2].coefficients[:, m+basis_order+1]]
    
    # How is this so small? The forcing was generated randomly?? 
    norm(A * x - b)

    # this should be very close to x as this is exactly how we solved for x  
    x2 = A \ b
    norm(x2 - x)


