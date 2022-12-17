function boundarycondition_modes(ω::AbstractFloat, bc::TractionBoundary, bearing::Bearing{2}, basis_order::Int)
    r = (bc.inner == true) ? bearing.r1 : bearing.r2
    return [
        pressure_traction_modes(ω, r, bearing.medium, basis_order);
        shear_traction_modes(ω, r, bearing.medium, basis_order)
    ]
end

function boundarycondition_modes(ω::AbstractFloat, bc::DisplacementBoundary, bearing::Bearing{2}, basis_order::Int)
    r = (bc.inner == true) ? bearing.r1 : bearing.r2
    return [
        pressure_displacement_modes(ω, r, bearing.medium, basis_order);
        shear_displacement_modes(ω, r, bearing.medium, basis_order)
    ]
end

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

function boundarycondition_systems(ω::AbstractFloat, bearing::Bearing{2}, bcs::AbstractBoundaryConditions, basis_order::Int)

    p1_j, p1_h, s1_j, s1_h = boundarycondition_modes(ω, bcs[1], bearing, basis_order);
    p2_j, p2_h, s2_j, s2_h = boundarycondition_modes(ω, bcs[2], bearing, basis_order);

    return M = map(1:(2basis_order+1)) do n
        vcat(
            hcat(p1_j[:,n],p1_h[:,n],s1_j[:,n],s1_h[:,n]),
            hcat(p2_j[:,n],p2_h[:,n],s2_j[:,n],s2_h[:,n])
        );
    end
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
