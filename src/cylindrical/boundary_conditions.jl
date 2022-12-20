function boundarycondition_mode(ω::AbstractFloat, bc::TractionBoundary, bearing::RollerBearing, basis_order::Int)
    r = (bc.inner == true) ? bearing.inner_radius : bearing.outer_radius
    return hcat(
        pressure_traction_mode(ω, r, bearing.medium, basis_order),
        shear_traction_mode(ω, r, bearing.medium, basis_order)
    )
end

function boundarycondition_mode(ω::AbstractFloat, bc::DisplacementBoundary, bearing::RollerBearing, basis_order::Int)
    r = (bc.inner == true) ? bearing.inner_radius : bearing.outer_radius
    return hcat(
        pressure_displacement_mode(ω, r, bearing.medium, basis_order),
        shear_displacement_mode(ω, r, bearing.medium, basis_order)
    )
end

function boundarycondition_system(ω::AbstractFloat, bearing::RollerBearing, bc1::AbstractBoundaryCondition, bc2::AbstractBoundaryCondition, basis_order::Int)

    first_boundary = boundarycondition_mode(ω, bc1, bearing, basis_order)
    second_boundary = boundarycondition_mode(ω, bc2, bearing, basis_order)

    return vcat(first_boundary, second_boundary)
end

# function ElasticWave(ω::T, bearing::RollerBearing, bcs::AbstractBoundaryConditions;
function ElasticWave(sim::BearingSimulation)

    ω = sim.ω
    bearing = sim.bearing

    kP = ω / bearing.medium.cp;
    kS = ω / bearing.medium.cs

    basis_order = sim.basis_order

    coefficients = map(-basis_order:basis_order) do m
        A = boundarycondition_system(ω, bearing, sim.boundarydata1.boundarytype, sim.boundarydata2.boundarytype, m)
        b = [
             sim.boundarydata1.fourier_modes[m+basis_order+1,:];
             sim.boundarydata2.fourier_modes[m+basis_order+1,:]
        ]

        x = A \ b

        relative_error = norm(A*x - b) / norm(b)
        if relative_error > sim.tol
            @warn "The relative error for the boundary conditions was $(relative_error) for (ω,basis_order) = $((ω,m))"
        end

        x
    end

    coefficients = transpose(hcat(coefficients...)) |> collect

    pressure_coefficients = coefficients[:,1:2] |> transpose
    shear_coefficients = coefficients[:,3:4] |> transpose

    φ = HelmholtzPotential{2}(kP,pressure_coefficients)
    ψ = HelmholtzPotential{2}(kS,shear_coefficients)

    return ElasticWave{2}(ω, bearing.medium, φ, ψ)
end
