function boundarycondition_pmode(ω::T, order::Int, bc::StressBoundary, bearing::Bearing{T,2}) where {T<:Number}

    r = (bc.inner == true) ? bearing.r1 : bearing.r2

    ρ = bearing.medium.ρ
    cp = bearing.medium.cp
    cs = bearing.medium.cs
    n = order
    ρs = ρ * cs^2 / r^2
    kP = ω / cp
    kS = ω / cs

    stress_column_p(J) = [
        -ρs * (2kP * r * J(n-1, kP*r) + (- 2n - 2n^2 + kS^2 * r^2) * J(n, kP*r)),
        im * n * ρs * (kP*r*J(n-1, kP*r) - 2J(n, kP*r) - kP * r * J(1 + n, kP*r))
    ]

    return hcat(stress_column_p(besselj), stress_column_p(hankelh1))
end

function boundarycondition_smode(ω::T, order::Int, bc::StressBoundary, bearing::Bearing{T,2}) where {T<:Number}

    r = (bc.inner == true) ? bearing.r1 : bearing.r2

    ρ = bearing.medium.ρ
    cp = bearing.medium.cp
    cs = bearing.medium.cs
    n = order
    ρs = ρ * cs^2 / r^2
    kS = ω / cs

    stress_column_s(J)=[
        im * n * ρs * (kS*r*J(n-1, kS*r) - 2J(n, kS*r) - kS*r*J(n+1, kS*r)),
        ρs * (2kS*r*J(n-1, kS*r) + (-2*n - 2*n^2 + kS^2*r^2) * J(n, kS*r))
    ]
    return hcat(stress_column_s(besselj), stress_column_s(hankelh1))
end

function boundarycondition_pmode(ω::T, order::Int, bc::DisplacementBoundary, bearing::Bearing{T,2}) where T<:Number

    r = (bc.inner == true) ? bearing.r1 : bearing.r2

    n = order
    cp = bearing.medium.cp
    kP = ω / cp

    displacement_column_p(J) = [
        kP/2 * (J(n-1, kP*r) - J(n+1, kP*r)), im * n * J(n, kP*r) / r
    ]

    return hcat(
        displacement_column_p(besselj), displacement_column_p(hankelh1)
    )
end

function boundarycondition_smode(ω::T, order::Int, bc::DisplacementBoundary, bearing::Bearing{T,2}) where {T<:Number}

    r = (bc.inner == true) ? bearing.r1 : bearing.r2

    n = order
    cs = bearing.medium.cs
    kS = ω / cs

    displacement_column_s(J) = [
        im * n * J(n, kS*r) / r, -kS/2 * (J(n-1, kS*r) - J(n+1, kS*r))
    ]
    return hcat(
        displacement_column_s(besselj), displacement_column_s(hankelh1)
    )
end

function system_matrix(ω::T, order::Int, bearing::Bearing{T,2}, bcs::AbstractBoundaryConditions) where T<:Number

    first_boundary = hcat(
        boundarycondition_pmode(ω, order, bcs[1], bearing), boundarycondition_smode(ω, order, bcs[1], bearing)
    )

    second_boundary = hcat(
        boundarycondition_pmode(ω, order, bcs[2], bearing), boundarycondition_smode(ω, order, bcs[2], bearing)
    )

    return vcat(first_boundary, second_boundary)
end

function mode_coefficients(ω::T, bearing::Bearing{T,2}, bcs::AbstractBoundaryConditions, forcing::Matrix{Complex{T}}) where T

    basis_order = Int( (size(forcing,1) - 1) / 2)

    coefficients = map(-basis_order:basis_order) do m
        A = system_matrix(ω, m, bearing, bcs)
        b = forcing[m+basis_order+1,:]
        A \ b
    end

    coefficients = transpose(hcat(coefficients...)) |> collect

    return coefficients
end
