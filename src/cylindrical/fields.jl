function displacement_modes(r::T, wave::ElasticWave{2}) where T <: AbstractFloat

    basis_order = wave.pressure.basis_order
    pmodes = hcat([
        pressure_displacement_mode(wave.ω, r, wave.medium, m) * wave.pressure.coefficients[:,m + basis_order + 1]
    for m = -basis_order:basis_order]...)

    basis_order = wave.shear.basis_order
    smodes = hcat([
        shear_displacement_mode(wave.ω, r, wave.medium, m) * wave.shear.coefficients[:,m + basis_order + 1]
    for m = -basis_order:basis_order]...)

    return transpose(pmodes + smodes) |> collect
end

function traction_modes(r::T, wave::ElasticWave{2}) where T <: AbstractFloat

    basis_order = wave.pressure.basis_order
    pmodes = hcat([
        pressure_traction_mode(wave.ω, r, wave.medium, m) * wave.pressure.coefficients[:,m + basis_order + 1]
    for m = -basis_order:basis_order]...)

    basis_order = wave.shear.basis_order
    smodes = hcat([
        shear_traction_mode(wave.ω, r, wave.medium, m) * wave.shear.coefficients[:,m + basis_order + 1]
    for m = -basis_order:basis_order]...)

    return transpose(pmodes + smodes) |> collect
end

function displacement(x::Vector{T}, wave::ElasticWave{2}) where T <: AbstractFloat

    r, θ = cartesian_to_radial_coordinates(x)
    # exps = exp.(im * θ .* (-basis_order:basis_order))

    basis_order = wave.pressure.basis_order
    modes = hcat([
        pressure_displacement_mode(wave.ω, r, wave.medium, m) .*  exp(im * θ * m)
    for m = -basis_order:basis_order]...)

    displace_p = modes * wave.pressure.coefficients[:]

    # displace_p = exp(im * m * θ) .* mode * wave.pressure.coefficients[m + basis_order + 1,:]

    basis_order = wave.shear.basis_order
    modes = hcat([
        shear_displacement_mode(wave.ω, r, wave.medium, m) .*  exp(im * θ * m)
    for m = -basis_order:basis_order]...)

    displace_s = modes * wave.shear.coefficients[:]

    return displace_p + displace_s
end

function traction(x::Vector{T}, wave::ElasticWave{2}) where T <: AbstractFloat

    r, θ = cartesian_to_radial_coordinates(x)

    basis_order = wave.pressure.basis_order
    modes = hcat([
        pressure_traction_mode(wave.ω, r, wave.medium, m) .*  exp(im * θ * m)
    for m = -basis_order:basis_order]...)

    traction_p = modes * wave.pressure.coefficients[:]

    basis_order = wave.shear.basis_order
    modes = hcat([
        shear_traction_mode(wave.ω, r, wave.medium, m) .*  exp(im * θ * m)
    for m = -basis_order:basis_order]...)

    traction_s = modes * wave.shear.coefficients[:]

    return traction_p + traction_s
end

# function pressure_mode(ω::AbstractFloat, bc::TractionBoundary, bearing::RollerBearing, basis_order::Int)
#     r = (bc.inner == true) ? bearing.r1 : bearing.r2
#     pressure_displacement_mode(ω, r, bearing.medium, basis_order)
# end


function pressure_displacement_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elasticity{2}, basis_order::Int)

    kP = ω / medium.cp;
    n = basis_order;

    bessel_modes(J::Function) = [kP/2 * (J(n-1, kP*r) - J(n+1, kP*r)), im * n * J(n, kP*r) / r]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function shear_displacement_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elasticity{2}, basis_order::Int)

    n = basis_order;
    cs = medium.cs
    kS = ω / cs

    bessel_modes(J::Function) = [im * n * J(n, kS*r) / r, -kS/2 * (J(n-1, kS*r) - J(n+1, kS*r))]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function pressure_traction_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elasticity{2}, basis_order::Int)

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

function shear_traction_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elasticity{2}, basis_order::Int)

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
