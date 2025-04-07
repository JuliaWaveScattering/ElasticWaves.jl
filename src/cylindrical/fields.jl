"""
    field(potential::HelmholtzPotential, x::AbstractVector)

Returns the value of the potentical at the point `x`
"""
function field(potential::HelmholtzPotential, x::AbstractVector{T}) where T

    k = potential.wavenumber
    r, θ = cartesian_to_radial_coordinates(x)

    coefs = permutedims(potential.coefficients, (2,1))
    exps = exp.(im * θ .* potential.modes)

    j = sum(coefs[:,1] .* besselj.(potential.modes,k*r) .* exps)
    h = sum(coefs[:,2] .* hankelh1.(potential.modes,k*r) .* exps)

   return j + h
end

function field(potential::HelmholtzPotential, bearing::RollerBearing; kws...)

    inner_circle = Circle(bearing.inner_radius)
    outer_circle = Circle(bearing.outer_radius)

    return field(potential, outer_circle;
        exclude_region = inner_circle,
        kws...
    )
end

function field_modes(wave::ElasticWave{2}, r::T, FT::FieldType) where T

    pmodes = hcat([
        pressure_field_mode(wave.ω, r, wave.medium, wave.potentials[1].modes[i], FT) * wave.potentials[1].coefficients[:,i]
    for i = eachindex(wave.potentials[1].modes)]...)

    smodes = hcat([
        shear_field_mode(wave.ω, r, wave.medium, wave.potentials[2].modes[i], FT) * wave.potentials[2].coefficients[:,i]
    for i = eachindex(wave.potentials[2].modes)]...)

    return transpose(pmodes + smodes) |> collect
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

    return  FrequencySimulationResult(reshape(field_mat, :, 1), x_vec, [wave.ω])
end


function field(wave::ElasticWave{2}, x::AbstractVector{T}, field_type::FieldType) where T 

    r, θ = cartesian_to_radial_coordinates(x)
    # exps = exp.(im * θ .* modes)

    modes_vec = [
        pressure_field_mode(wave.ω, r, wave.medium, m, field_type) .*  exp(im * θ * m)
    for m = wave.potentials[1].modes]

    modes_matrix = hcat(modes_vec...)

    field_p = modes_matrix * wave.potentials[1].coefficients[:]

    # displace_p = exp(im * m * θ) .* mode * wave.potentials[1].coefficients[m + basis_order + 1,:]

    modes_matrix = hcat([
        shear_field_mode(wave.ω, r, wave.medium, m, field_type) .*  exp(im * θ * m)
    for m = wave.potentials[2].modes]...)

    field_s = modes_matrix * wave.potentials[2].coefficients[:]

    return field_p + field_s
end


function pressure_field_mode(ω::T, r::Union{T,Complex{T}}, medium::Elastic{2}, mode::Int, ::DisplacementType) where T

    kP = ω / medium.cp;
    n = mode;

    bessel_modes(J::Function) = [kP/2 * (J(n-1, kP*r) - J(n+1, kP*r)), im * n * J(n, kP*r) / r]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function shear_field_mode(ω::T, r::Union{T,Complex{T}}, medium::Elastic{2}, mode::Int, ::DisplacementType) where T

    n = mode;
    cs = medium.cs
    kS = ω / cs

    bessel_modes(J::Function) = [im * n * J(n, kS*r) / r, -kS/2 * (J(n-1, kS*r) - J(n+1, kS*r))]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function pressure_field_mode(ω::T, r::Union{T,Complex{T}}, medium::Elastic{2}, mode::Int, ::TractionType) where T

    ρ = medium.ρ
    n = mode;
    cp = medium.cp; cs = medium.cs
    kP = ω / cp; kS = ω / cs

    ρs = ρ * cs^2 / r^2
    bessel_modes(J::Function) = [
        -ρs * (2kP * r * J(n-1, kP*r) + (- 2n - 2n^2 + kS^2 * r^2) * J(n, kP*r)),
        im * n * ρs * (kP*r*J(n-1, kP*r) - 2J(n, kP*r) - kP * r * J(1 + n, kP*r))
    ]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function shear_field_mode(ω::T, r::Union{T,Complex{T}}, medium::Elastic{2}, mode::Int, ::TractionType) where T

    ρ = medium.ρ
    n = mode;
    cp = medium.cp; cs = medium.cs
    kP = ω / cp; kS = ω / cs

    ρs = ρ * cs^2 / r^2
    bessel_modes(J::Function) = [
        im * n * ρs * (kS*r*J(n-1, kS*r) - 2J(n, kS*r) - kS*r*J(n+1, kS*r)),
        ρs * (2kS*r*J(n-1, kS*r) + (-2*n - 2*n^2 + kS^2 * r^2) * J(n, kS*r))
    ]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end
