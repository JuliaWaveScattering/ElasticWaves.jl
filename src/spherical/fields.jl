# function field(wave::ElasticWave{3}, x::AbstractVector{T}, field_type::FieldType) where T <: AbstractFloat

#     # pressure_regular_basis, shearΦ_regular_basis, shearχ_regular_basis

#     ω = 1.1
#     basis_order = 8

#     # medium = Elastic(3; ρ = 0.5, cp = 1.1, cs = 0.9)


# end

## TractionType fields

function pressure_regular_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::TractionType) where T
    
    rθφ = cartesian_to_radial_coordinates(x)
    r, θ, φ = rθφ

    μ = medium.cs^2 * medium.ρ
    λ = medium.cp^2 * medium.ρ - 2μ

    kp = ω / medium.cp;
    kpr = kp * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [sbesselj(l, kpr) for l = 0:basis_order]
    djs = [diffsbesselj(l, kpr) for l = 0:basis_order]
    ddjs  = [diff2sbesselj(l, ksr) for l = 0:basis_order]

    lm_to_n = lm_to_spherical_harmonic_index
    
    # need to transform the 3D vectors from a spherical to a Cartesian coordinate basis
    M = spherical_to_cartesian_transform(rθφ)

    ps = [
        transpose(
            M * [
                kp^2 * Ys[lm_to_n(l,m)] * (-js[l+1]* λ + 2ddjs[l+1]* μ),
                (2dYs[lm_to_n(l,m)] * (-js[l+1]+ djs[l+1] * kpr) * μ) / r^2,
                (2im * cscθ * m * (-js[l+1] + djs[l+1] * kpr) * Ys[lm_to_n(l,m)] * μ) / r^2
            ]
        )
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end


function shearΦ_regular_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::TractionType) where T
    
    rθφ  = cartesian_to_radial_coordinates(x)
    r, θ, φ  = rθφ

    μ = medium.cs^2 * medium.ρ

    ks = ω / medium.cs;
    ksr = ks * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [sbesselj(l, ksr) for l = 0:(basis_order+1)]
    djs = [diffsbesselj(l, ksr) for l = 0:basis_order]
    ddjs  = [diff2sbesselj(l, ksr) for l = 0:basis_order]
    # ddjs = [
    #     2js[l+2] / ksr + (l^2 - l - ksr^2) * js[l+1] / ksr^2
    # for l = 0:basis_order]

    lm_to_n = lm_to_spherical_harmonic_index

    # need to transform the 3D vectors from a spherical to a Cartesian coordinate basis
    M = spherical_to_cartesian_transform(rθφ)

    ps = [
        M * [
            2ks^2 * (3ddjs[l+1] + js[l+1] + (dddjs[l+1] + djs[l+1]) * ksr) * Ys[lm_to_n(l,m)] * μ, 
            (dYs[lm_to_n(l,m)] * ( 2ksr * (djs[l+1] + ddjs[l+1] * ksr) + js[l+1] * (-2 + ks^2 * r^2)) * μ) / r^2,
            (im * cscθ * m * (2ksr * (djs[l+1] + ddjs[l+1] * ksr) + js[l+1] * (-2 + ks^2 * r^2)) * Ys[lm_to_n(l,m)] * μ) / r^2
        ] |> transpose
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end


function shearχ_regular_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::TractionType) where T
    
    rθφ  = cartesian_to_radial_coordinates(x)
    r, θ, φ  = rθφ

    μ = medium.cs^2 * medium.ρ
    
    ks = ω / medium.cs;
    ksr = ks * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [sbesselj(l, ksr) for l = 0:(basis_order+1)]

    lm_to_n = lm_to_spherical_harmonic_index

    # the components of the vectors below are in a spherical coordinate basis. That is the vectors below have no radial component
    # need to transform the 3D vectors from a spherical to a Cartesian coordinate basis
    M = spherical_to_cartesian_transform(rθφ)

    ps = [
        M * [
            zero(Complex{T}), 
            (im * cscθ * ks * m * (-js[l+1] + djs[l+1] * ksr) * Ys[lm_to_n(l,m)] * μ) / r,
            (dYs[lm_to_n(l,m)] * ks * (js[l+1] - djs[l+1] * ksr) * μ) / r
        ] |> transpose
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end


## DisplacementType fields

function pressure_regular_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::DisplacementType) where T
    
    rθφ = cartesian_to_radial_coordinates(x)
    r, θ, φ = rθφ

    kp = ω / medium.cp;
    kpr = kp * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [sbesselj(l, kpr) for l = 0:basis_order]
    djs = [diffsbesselj(l, kpr) for l = 0:basis_order]

    lm_to_n = lm_to_spherical_harmonic_index
    
    # need to transform the 3D vectors from a spherical to a Cartesian coordinate basis
    M = spherical_to_cartesian_transform(rθφ)

    ps = [
        transpose(
            M * [
                kp * Ys[lm_to_n(l,m)] * djs[l+1], 
                js[l+1] * dYs[lm_to_n(l,m)] / r, 
                im * m * cscθ * js[l+1] * Ys[lm_to_n(l,m)] / r
            ]
        )        
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end    

function pressure_outgoing_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::DisplacementType) where T
    
    rθφ = cartesian_to_radial_coordinates(x)
    r, θ, φ = rθφ

    kp = ω / medium.cp;
    kpr = kp * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [shankelh1(l, kpr) for l = 0:basis_order]
    djs = [diffshankelh1(l, kpr) for l = 0:basis_order]

    lm_to_n = lm_to_spherical_harmonic_index
    
    # need to transform the 3D vectors from a spherical to a Cartesian coordinate basis
    M = spherical_to_cartesian_transform(rθφ)

    ps = [
        transpose(
            M * [
                kp * Ys[lm_to_n(l,m)] * djs[l+1], 
                js[l+1] * dYs[lm_to_n(l,m)] / r, 
                im * m * cscθ * js[l+1] * Ys[lm_to_n(l,m)] / r
            ]
        )        
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end    

function shearΦ_regular_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::DisplacementType) where T
    
    rθφ  = cartesian_to_radial_coordinates(x)
    r, θ, φ  = rθφ

    ks = ω / medium.cs;
    ksr = ks * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [sbesselj(l, ksr) for l = 0:(basis_order+1)]
    djs = [diffsbesselj(l, ksr) for l = 0:basis_order]
    ddjs  = [diff2sbesselj(l, ksr) for l = 0:basis_order]
    # ddjs = [
    #     2js[l+2] / ksr + (l^2 - l - ksr^2) * js[l+1] / ksr^2
    # for l = 0:basis_order]

    lm_to_n = lm_to_spherical_harmonic_index

    # need to transform the 3D vectors from a spherical to a Cartesian coordinate basis
    M = spherical_to_cartesian_transform(rθφ)

    ps = [
        M * [
            ks * Ys[lm_to_n(l,m)] * (ksr * js[l+1] + 2djs[l+1] + ksr * ddjs[l+1]), 
            (js[l+1] + ksr * djs[l+1]) * dYs[lm_to_n(l,m)] / r,
            im*m * cscθ * (js[l+1] + ksr * djs[l+1]) * Ys[lm_to_n(l,m)] / r
        ] |> transpose
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end    

function shearΦ_outgoing_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::DisplacementType) where T
    
    rθφ  = cartesian_to_radial_coordinates(x)
    r, θ, φ  = rθφ

    ks = ω / medium.cs;
    ksr = ks * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [shankelh1(l, ksr) for l = 0:(basis_order+1)]
    djs = [diffshankelh1(l, ksr) for l = 0:basis_order]
    ddjs  = [diff2shankelh1(l, ksr) for l = 0:basis_order]

    lm_to_n = lm_to_spherical_harmonic_index

    # need to transform the 3D vectors from a spherical to a Cartesian coordinate basis
    M = spherical_to_cartesian_transform(rθφ)

    ps = [
        M * [
            ks * Ys[lm_to_n(l,m)] * (ksr * js[l+1] + 2djs[l+1] + ksr * ddjs[l+1]), 
            (js[l+1] + ksr * djs[l+1]) * dYs[lm_to_n(l,m)] / r,
            im*m * cscθ * (js[l+1] + ksr * djs[l+1]) * Ys[lm_to_n(l,m)] / r
        ] |> transpose
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end    


function shearχ_regular_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::DisplacementType) where T
    
    rθφ  = cartesian_to_radial_coordinates(x)
    r, θ, φ  = rθφ
    
    ks = ω / medium.cs;
    ksr = ks * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [sbesselj(l, ksr) for l = 0:(basis_order+1)]

    lm_to_n = lm_to_spherical_harmonic_index

    # the components of the vectors below are in a spherical coordinate basis. That is the vectors below have no radial component
    # need to transform the 3D vectors from a spherical to a Cartesian coordinate basis
    M = spherical_to_cartesian_transform(rθφ)

    ps = [
        M * [
            zero(Complex{T}), 
            im*m * ks * cscθ * js[l+1] * Ys[lm_to_n(l,m)],
            - ks * js[l+1] * dYs[lm_to_n(l,m)]
        ] |> transpose
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end    

function shearχ_outgoing_basis(ω::AbstractFloat, x::AbstractVector{T}, medium::Elastic{3}, basis_order::Int, ::DisplacementType) where T
    
    rθφ  = cartesian_to_radial_coordinates(x)
    r, θ, φ  = rθφ
    
    ks = ω / medium.cs;
    ksr = ks * r
    cscθ = csc(θ)

    Ys = spherical_harmonics(basis_order, θ, φ)
    dYs = spherical_harmonics_dθ(basis_order, θ, φ)
    
    js = [shankelh1(l, ksr) for l = 0:(basis_order+1)]

    lm_to_n = lm_to_spherical_harmonic_index

    # the components of the vectors below are in a spherical coordinate basis. That is the vectors below have no radial component
    # need to transform the 3D vectors from a spherical to a Cartesian coordinate basis
    M = spherical_to_cartesian_transform(rθφ)

    ps = [
        M * [
            zero(Complex{T}), 
            im*m * ks * cscθ * js[l+1] * Ys[lm_to_n(l,m)],
            - ks * js[l+1] * dYs[lm_to_n(l,m)]
        ] |> transpose
    for l = 0:basis_order for m = -l:l]

    return vcat(ps...)
end

