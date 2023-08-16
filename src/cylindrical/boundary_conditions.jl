function boundarycondition_mode(ω::AbstractFloat, bc::BoundaryCondition, bearing::RollerBearing, basis_order::Int)
    r = (bc.inner == true) ? bearing.inner_radius : bearing.outer_radius
    return hcat(
        pressure_field_mode(ω, r, bearing.medium, basis_order, bc.fieldtype),
        shear_field_mode(ω, r, bearing.medium, basis_order, bc.fieldtype)
    )
end


function boundarycondition_system(ω::AbstractFloat, bearing::RollerBearing, bc1::BoundaryCondition, bc2::BoundaryCondition, basis_order::Int)

    first_boundary = boundarycondition_mode(ω, bc1, bearing, basis_order)
    second_boundary = boundarycondition_mode(ω, bc2, bearing, basis_order)

    return vcat(first_boundary, second_boundary)
end

function ElasticWave(sim::BearingSimulation)

    ω = sim.ω
    bearing = sim.bearing

    basis = sim.boundarybasis
    T = typeof(ω)

    kP = ω / bearing.medium.cp;
    kS = ω / bearing.medium.cs

    basis_order = sim.basis_order

    coefficients = [
        zeros(Complex{T}, 4)
    for m = -basis_order:basis_order]

    mode_errors = zeros(T, 2basis_order + 1)

     for m in -basis_order:basis_order
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

        coefficients[m + basis_order + 1] = x
        mode_errors[m + basis_order + 1] = relative_error
    end

    coefficients = transpose(hcat(coefficients...)) |> collect

    pressure_coefficients = coefficients[:,1:2] |> transpose
    shear_coefficients = coefficients[:,3:4] |> transpose

    φ = HelmholtzPotential{2}(bearing.medium.cp, kP, pressure_coefficients)
    ψ = HelmholtzPotential{2}(bearing.medium.cs, kS, shear_coefficients)

    return ElasticWave(ω, bearing.medium, φ, ψ; mode_errors = mode_errors)
end
