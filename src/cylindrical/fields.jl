function field_modes(wave::ElasticWave{2}, r::T, FT::FieldType) where T <: AbstractFloat

    basis_order = wave.pressure.basis_order
    pmodes = hcat([
        pressure_field_mode(wave.ω, r, wave.medium, m, FT) * wave.pressure.coefficients[:,m + basis_order + 1]
    for m = -basis_order:basis_order]...)

    basis_order = wave.shear.basis_order
    smodes = hcat([
        shear_field_mode(wave.ω, r, wave.medium, m, FT) * wave.shear.coefficients[:,m + basis_order + 1]
    for m = -basis_order:basis_order]...)

    return transpose(pmodes + smodes) |> collect
end


function displacement(wave::ElasticWave, x::Vector{T}) where T <: AbstractFloat
    field(wave, x, DisplacementType())
end

function traction(wave::ElasticWave, x::Vector{T}) where T <: AbstractFloat
    field(wave, x, TractionType())
end

function field(wave::ElasticWave{2}, bearing::RollerBearing, fieldtype::FieldType; kws...)

    inner_circle = Circle(bearing.inner_radius)
    outer_circle = Circle(bearing.outer_radius)

    return field(wave, outer_circle, fieldtype;
        exclude_region = inner_circle,
        kws...
    )
end

function field(wave::ElasticWave{2}, sh::Shape, fieldtype::FieldType; kws...)

    x_vec, inds = points_in_shape(sh; kws...)
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.

    xs = x_vec[inds]
    field_mat = [SVector(0.0+0.0im, 0.0+0.0im) for x in x_vec]

    fs = [field(wave,x,fieldtype) for x in xs];
    field_mat[inds] = fs

    return  FrequencySimulationResult(reshape(field_mat, :, 1), x_vec, [ω])
end

function field(wave::ElasticWave{2}, x::AbstractVector{T}, FT::FieldType) where T <: AbstractFloat

    r, θ = cartesian_to_radial_coordinates(x)
    # exps = exp.(im * θ .* (-basis_order:basis_order))

    basis_order = wave.pressure.basis_order
    modes = hcat([
        pressure_field_mode(wave.ω, r, wave.medium, m, FT) .*  exp(im * θ * m)
    for m = -basis_order:basis_order]...)

    field_p = modes * wave.pressure.coefficients[:]

    # displace_p = exp(im * m * θ) .* mode * wave.pressure.coefficients[m + basis_order + 1,:]

    basis_order = wave.shear.basis_order
    modes = hcat([
        shear_field_mode(wave.ω, r, wave.medium, m, FT) .*  exp(im * θ * m)
    for m = -basis_order:basis_order]...)

    field_s = modes * wave.shear.coefficients[:]

    return field_p + field_s
end


function pressure_field_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elasticity{2}, basis_order::Int, ::DisplacementType)

    kP = ω / medium.cp;
    n = basis_order;

    bessel_modes(J::Function) = [kP/2 * (J(n-1, kP*r) - J(n+1, kP*r)), im * n * J(n, kP*r) / r]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function shear_field_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elasticity{2}, basis_order::Int, ::DisplacementType)

    n = basis_order;
    cs = medium.cs
    kS = ω / cs

    bessel_modes(J::Function) = [im * n * J(n, kS*r) / r, -kS/2 * (J(n-1, kS*r) - J(n+1, kS*r))]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function pressure_field_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elasticity{2}, basis_order::Int, ::TractionType)

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

function shear_field_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elasticity{2}, basis_order::Int, ::TractionType)

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
