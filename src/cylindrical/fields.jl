function displacement(x::Vector{V}, ElasticWave::ElasticWave{2}) where V <: AbstractVector

    r, θ = cartesian_to_radial_coordinates(x)

    pjs, phs = pressure_displacement_modes(ElasticWave.ω, r, ElasticWave.basis_order, ElasticWave.medium)
    sjs, shs = shear_displacement_modes(ElasticWave.ω, r, ElasticWave.basis_order, ElasticWave.medium)

    return modes
end

function pressure_displacement_modes(ω::AbstractFloat, r::AbstractFloat, basis_order::Int, medium::Elasticity{2})

    cp = medium.cp
    kP = ω / cp;

    bessel_modes(J::Function) = map(-basis_order:basis_order) do n
        kP/2 * (J(n-1, kP*r) - J(n+1, kP*r)), im * n * J(n, kP*r) / r
    end

    besselj_modes = hcat(bessel_modes(besselj)...)
    hankelh1_modes = hcat(bessel_modes(hankelh1)...)

    return [besselj_modes, hankelh1_modes]

end

function shear_displacement_modes(ω::AbstractFloat, r::AbstractFloat, basis_order::Int, medium::Elasticity{2})

    cs = medium.cs
    kS = ω / cs

    bessel_modes(J::Function) = map(-basis_order:basis_order) do n
        im * n * J(n, kS*r) / r, -kS/2 * (J(n-1, kS*r) - J(n+1, kS*r))
    end

    besselj_modes = hcat(bessel_modes(besselj)...)
    hankelh1_modes = hcat(bessel_modes(hankelh1)...)

    return [besselj_modes, hankelh1_modes]
end

function pressure_traction_modes(ω::AbstractFloat, r::AbstractFloat, basis_order::Int, medium::Elasticity{2})

    ρ = medium.ρ
    cp = medium.cp; cs = medium.cs

    kP = ω / cp; kS = ω / cs

    bessel_modes(J::Function) = map(-basis_order:basis_order) do n
        [-ρ * cs^2 / r^2 * (2kP * r * J(n-1, kP*r) + (- 2n - 2n^2 + kS^2 * r^2) * J(n, kP*r)),
        im * n * ρ * cs^2 / r^2 * (kP*r*J(n-1, kP*r) - 2J(n, kP*r) - kP * r * J(1 + n, kP*r))]
    end

    besselj_modes = hcat(bessel_modes(besselj)...)
    hankelh1_modes = hcat(bessel_modes(hankelh1)...)

    # the reshape below means that M[:,:,n] * [a[n], b[n]] returns the stress for basis_order n
    # M = reshape(vcat(besselmodes, hankelmodes), (2, 2,:))
    return [besselj_modes, hankelh1_modes]

end

function shear_traction_modes(ω::AbstractFloat, r::AbstractFloat, basis_order::Int, medium::Elasticity{2})

    ρ = medium.ρ
    cp = medium.cp; cs = medium.cs

    kP = ω / cp; kS = ω / cs

    bessel_modes(J::Function) = map(-basis_order:basis_order) do n
        [im * n * ρ * cs^2 / r^2 * (kS*r*J(n-1, kS*r) - 2J(n, kS*r) - kS*r*J(n+1, kS*r)),
        ρ * cs^2 / r^2 * (2kS*r*J(n-1, kS*r) + (-2*n - 2*n^2 + kS^2 * r^2) * J(n, kS*r))]
    end

    besselj_modes = hcat(bessel_modes(besselj)...)
    hankelh1_modes = hcat(bessel_modes(hankelh1)...)

    # vcat(besselmodes, hankelmodes)
    return [besselj_modes, hankelh1_modes]
end

function pressure_displacement_mode(ω::AbstractFloat, r::AbstractFloat, basis_order::Int, medium::Elasticity{2})

    kP = ω / medium.cp;
    n = basis_order;

    bessel_modes(J::Function) = [kP/2 * (J(n-1, kP*r) - J(n+1, kP*r)), im * n * J(n, kP*r) / r]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function shear_displacement_mode(ω::AbstractFloat, r::AbstractFloat, basis_order::Int, medium::Elasticity{2})

    n = basis_order;
    cs = medium.cs
    kS = ω / cs

    bessel_modes(J::Function) = [im * n * J(n, kS*r) / r, -kS/2 * (J(n-1, kS*r) - J(n+1, kS*r))]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function pressure_traction_mode(ω::AbstractFloat, r::AbstractFloat, basis_order::Int, medium::Elasticity{2})

    ρ = medium.ρ
    n = basis_order;
    cp = medium.cp; cs = medium.cs
    kP = ω / cp; kS = ω / cs

    ρs = ρ * cs^2 / r^2
    bessel_modes(J::Function) = [
        -ρs * (2kP * r * J(n-1, kP*r) + (- 2n - 2n^2 + kS^2 * r^2) * J(n, kP*r)),
        im * n * ρs * (kP*r*J(n-1, kP*r) - 2J(n, kP*r) - kP * r * J(1 + n, kP*r))
    ]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function shear_traction_mode(ω::AbstractFloat, r::AbstractFloat, basis_order::Int, medium::Elasticity{2})

    ρ = medium.ρ
    n = basis_order;
    cp = medium.cp; cs = medium.cs
    kP = ω / cp; kS = ω / cs

    ρs = ρ * cs^2 / r^2
    bessel_modes(J::Function) = [
        im * n * ρs * (kS*r*J(n-1, kS*r) - 2J(n, kS*r) - kS*r*J(n+1, kS*r)),
        ρs * (2kS*r*J(n-1, kS*r) + (-2*n - 2*n^2 + kS^2 * r^2) * J(n, kS*r))
    ]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end



function pressure_potential(basis_order::Int; ω::T,r::T,θ::T,bearing::Bearing{2},fp_coefficients::Vector,fs_coefficients::Vector) where {T<:Number}
    n = basis_order
    f = coes(n;ω,bearing,fp_coefficients, fs_coefficients)
    kp = ω/bearing.medium.cp
    return (f[1]*besselj(n,kp*r) + f[2]*hankelh1(n,kp*r))*ℯ^(im*n*θ)
end

function shear_potential(basis_order::Int; ω::T,r::T,θ::T,bearing::Bearing{2},fp_coefficients::Vector,fs_coefficients::Vector) where {T<:Number}
    n = basis_order
    f = coes(n;ω,bearing,fp_coefficients, fs_coefficients)
    ks = ω/bearing.medium.cs
    return (f[3]*besselj(n,ks*r) + f[4]*hankelh1(n,ks*r))*ℯ^(im*n*θ)
end

function displacement_mode(basis_order::Int;ω::T,r::T,θ::T,bearing::Bearing{2},fp_coefficients::Vector,fs_coefficients::Vector) where {T<:Number}
    n = basis_order
    kp = ω/bearing.medium.cp
    ks = ω/bearing.medium.cs
    f = coes(n;ω,bearing,fp_coefficients,fs_coefficients)

    u1 = cos(θ) * (
            kp * (f[1]*diffbesselj(n,kp*r) + f[2]*diffhankelh1(n,kp*r)) +
            ((im*n)/r) * potentialψ(n;ω,r,θ,bearing,fp_coefficients,fs_coefficients)
        ) - sin(θ) * (
            ((im*n)/r) * potentialϕ(n;ω,r,θ,bearing,fp_coefficients,fs_coefficients) -ks*(f[3]*diffbesselj(n,ks*r)+f[4]*diffhankelh1(n,ks*r))
        )

    u2 = sin(θ) * (
            kp * (f[1]*diffbesselj(n,kp*r)+f[2]*diffhankelh1(n,kp*r)) +
            ((im*n)/r) * potentialψ(n;ω,r,θ,bearing,fp_coefficients,fs_coefficients)
        ) + cos(θ) * (
            ((im*n)/r) * potentialϕ(n;ω,r,θ,bearing,fp_coefficients,fs_coefficients) -ks*(f[3]*diffbesselj(n,ks*r)+f[4]*diffhankelh1(n,ks*r))
        )

    return [u1;u2]*ℯ^(im*n*θ)
end

function displacement(basis_order::Int;ω::T,r::T,θ::T,bearing::Bearing{2},fp_coefficients::Vector,fs_coefficients::Vector) where {T<:Number}

    u=zeros(2,1)

    for i in 1:basis_order
        u=u+displacement_mode(convert(Int,i-(basis_order+1)/2);ω,r,θ,bearing,fp_coefficients, fs_coefficients)
    end

    return u
end
