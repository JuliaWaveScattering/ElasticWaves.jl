function stress_matrix_pmode(ω::T, order::Int, r::T, bearing::Bearing{T,2}) where {T<:Number}
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

function stress_matrix_smode(ω::T, order::Int, r::T, bearing::Bearing{T,2}) where {T<:Number}
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

function displacement_matrix_pmode(ω::T, order::Int, r::T, bearing::Bearing{T,2}) where T<:Number
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

function displacement_matrix_smode(ω::T, order::Int, r::T, bearing::Bearing{T,2}) where {T<:Number}
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

function stress_matrix_half(ω::T,  order::Int, r::T, bearing::Bearing{T,2}) where T<:Number
    return hcat(
        stress_matrix_pmode(ω,order,r,bearing), stress_matrix_smode(ω,order,r,bearing)
    )
end

function displacement_matrix_half(ω::T,  order::Int, r::T, bearing::Bearing{T,2}) where T<:Number
    return hcat(
        displacment_matrix_pmode(ω,order,r,bearing), displacment_matrix_smode(ω,order,r,bearing)
    )
end

function stress_matrix_full(order::Int; ω::T, bearing::Bearing{T,2}) where T<:Number
    return vcat(
        stress_matrix_half(ω,  order, bearing.r1, bearing),
        stress_matrix_half(ω,  order, bearing.r2, bearing)
    )
end


function forcings(order::Int;fp_coefficients::Vector, fs_coefficients::Vector)
    fourier_coes=transpose(hcat(fp_coefficients,fs_coefficients))
    n = order
    return fourier_coes[:,convert(Int,n+(1+size(fourier_coes,2))/2)]
end

function coes(order::Int; ω::T, bearing::Bearing{T,2}, fp_coefficients::Vector,fs_coefficients::Vector) where {T<:Number}
    n=order
    f=forcings(n; fp_coefficients, fs_coefficients)
    f_new = [f[1],f[2],0.0,0.0]


    return stress_matrix_full(n; ω, bearing) \ f_new
end
