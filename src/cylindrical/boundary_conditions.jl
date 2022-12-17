function boundarycondition_mode(ω::AbstractFloat, bc::TractionBoundary, bearing::Bearing{2}, basis_order::Int)
    r = (bc.inner == true) ? bearing.r1 : bearing.r2
    return hcat(
        pressure_traction_mode(ω, r, bearing.medium, basis_order),
        shear_traction_mode(ω, r, bearing.medium, basis_order)
    )
end

function boundarycondition_mode(ω::AbstractFloat, bc::DisplacementBoundary, bearing::Bearing{2}, basis_order::Int)
    r = (bc.inner == true) ? bearing.r1 : bearing.r2
    return hcat(
        pressure_displacement_mode(ω, r, bearing.medium, basis_order),
        shear_displacement_mode(ω, r, bearing.medium, basis_order)
    )
end

function boundarycondition_system(ω::AbstractFloat, bearing::Bearing{2}, bcs::AbstractBoundaryConditions, basis_order::Int)

    first_boundary = boundarycondition_mode(ω, bcs[1], bearing, basis_order)
    second_boundary = boundarycondition_mode(ω, bcs[2], bearing, basis_order)

    return vcat(first_boundary, second_boundary)
end

function ElasticWave(ω::T, bearing::Bearing{2}, bcs::AbstractBoundaryConditions, forcing_modes::Matrix{C};
        tol::T = eps(T)^(1/4)
    ) where {T, C <: Complex{T}}

    kP = ω / bearing.medium.cp;
    kS = ω / bearing.medium.cs

    basis_order = Int( (size(forcing_modes,1) - 1) / 2)

    coefficients = map(-basis_order:basis_order) do m
        A = boundarycondition_system(ω, bearing, bcs, m)
        b = forcing_modes[m+basis_order+1,:]
        x = A \ b

        relative_error = norm(A*x - b) / norm(b)
        if relative_error > tol
            @warn "The relative error for the boundary conditions was $(relative_error) for (ω,basis_order) = $((ω,m))"
        end

        x
    end

    coefficients = transpose(hcat(coefficients...)) |> collect

    pressure_coefficients = coefficients[:,1:2] |> transpose
    shear_coefficients = coefficients[:,3:4] |> transpose

    φ = HelmholtzPotential{2}(kP,pressure_coefficients)
    ψ = HelmholtzPotential{2}(kS,shear_coefficients)

    return ElasticWave{2}(ω,bearing.medium, φ, ψ)
end
