function field_modes(wave::ElasticWave{2}, r::T, FT::FieldType) where T <: AbstractFloat

    basis_order = wave.potentials[1].basis_order
    pmodes = hcat([
        pressure_field_mode(wave.ω, r, wave.medium, m, FT) * wave.potentials[1].coefficients[:,m + basis_order + 1]
    for m = -basis_order:basis_order]...)

    basis_order = wave.potentials[2].basis_order
    smodes = hcat([
        shear_field_mode(wave.ω, r, wave.medium, m, FT) * wave.potentials[2].coefficients[:,m + basis_order + 1]
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

    return  FrequencySimulationResult(reshape(field_mat, :, 1), x_vec, [wave.ω])
end

function field(wave::ElasticWave{3}, x::AbstractVector{T}, field_type::FieldType) where T <: AbstractFloat

    # pressure_field_basis, shearΦ_field_basis, shearχ_field_basis

    ω = 1.1
    basis_order = 8

    medium = Elastic(3; ρ = 0.5, cp = 1.1, cs = 0.9)


end



function pressure_field_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::DisplacementType) where T
    
    r, θ, φ  = cartesian_to_radial_coordinates(x)
    kp = ω / medium.cp;
    kpr = kP * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [sbesselj(l, kpr) for l = 0:basis_order]
    djs = [diffsbesselj(l, kpr) for l = 0:basis_order]

    lm_to_n = lm_to_spherical_harmonic_index
    
    # the components of the vectors below are in a spherical coordinate basis
    ps = [
        [kp * Ys[lm_to_n(l,m)] * djs[l+1], js[l+1] * dYs[lm_to_n(l,m)], im * m * cscθ * js[l+1] * Ys[lm_to_n(l,m)]] |> transpose
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end    

function shearΦ_field_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::DisplacementType) where T
    
    r, θ, φ  = cartesian_to_radial_coordinates(x)
    ks = ω / medium.cs;
    ksr = ks * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [sbesselj(l, ksr) for l = 0:(basis_order+1)]
    djs = [diffsbesselj(l, ksr) for l = 0:basis_order]
    ddjs = [
        2js[l+2] / ksr + (l^2 - l - ksr^2) * js[l+1] / ksr^2
    for l = 0:basis_order]

    lm_to_n = lm_to_spherical_harmonic_index

    # the components of the vectors below are in a spherical coordinate basis
    ps = [
        [
            ks * Ys[lm_to_n(l,m)] * (ksr * js[l+1] + 2djs[l+1] + ksr * ddjs), 
            (js[l+1] + ksr * djs[l+1]) * dYs[lm_to_n(l,m)] / r,
            im*m * cscθ * (js[l+1] + ksr * djs[l+1]) * Ys[lm_to_n(l,m)] / r
        ] |> transpose
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end    

function shearχ_field_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::DisplacementType) where T
    
    r, θ, φ  = cartesian_to_radial_coordinates(x)
    ks = ω / medium.cs;
    ksr = ks * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [sbesselj(l, ksr) for l = 0:(basis_order+1)]
    djs = [diffsbesselj(l, ksr) for l = 0:basis_order]
    ddjs = [
        2js[l+2] / ksr + (l^2 - l - ksr^2) * js[l+1] / ksr^2
    for l = 0:basis_order]

    lm_to_n = lm_to_spherical_harmonic_index

    # the components of the vectors below are in a spherical coordinate basis. That is the vectors below have no radial component
    ps = [
        [
            zero(Complex{T}), 
            im*m * ksr * cscθ * js[l+1] * Ys[lm_to_n(l,m)],
            - ksr * js[l+1] * dYs[lm_to_n(l,m)]
        ] |> transpose
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end    


function field(wave::ElasticWave{2}, x::AbstractVector{T}, field_type::FieldType) where T <: AbstractFloat

    r, θ = cartesian_to_radial_coordinates(x)
    # exps = exp.(im * θ .* (-basis_order:basis_order))

    basis_order = wave.potentials[1].basis_order
    modes_vec = [
        pressure_field_mode(wave.ω, r, wave.medium, m, field_type) .*  exp(im * θ * m)
    for m = -basis_order:basis_order]

    modes_matrix = hcat(modes_vec...)

    field_p = modes_matrix * wave.potentials[1].coefficients[:]

    # displace_p = exp(im * m * θ) .* mode * wave.potentials[1].coefficients[m + basis_order + 1,:]

    basis_order = wave.potentials[2].basis_order
    modes = hcat([
        shear_field_mode(wave.ω, r, wave.medium, m, field_type) .*  exp(im * θ * m)
    for m = -basis_order:basis_order]...)

    field_s = modes * wave.potentials[2].coefficients[:]

    return field_p + field_s
end


function pressure_field_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elastic{2}, basis_order::Int, ::DisplacementType)

    kP = ω / medium.cp;
    n = basis_order;

    bessel_modes(J::Function) = [kP/2 * (J(n-1, kP*r) - J(n+1, kP*r)), im * n * J(n, kP*r) / r]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function shear_field_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elastic{2}, basis_order::Int, ::DisplacementType)

    n = basis_order;
    cs = medium.cs
    kS = ω / cs

    bessel_modes(J::Function) = [im * n * J(n, kS*r) / r, -kS/2 * (J(n-1, kS*r) - J(n+1, kS*r))]

    return hcat(bessel_modes(besselj), bessel_modes(hankelh1))
end

function pressure_field_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elastic{2}, basis_order::Int, ::TractionType)

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

function shear_field_mode(ω::AbstractFloat, r::AbstractFloat, medium::Elastic{2}, basis_order::Int, ::TractionType)

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
