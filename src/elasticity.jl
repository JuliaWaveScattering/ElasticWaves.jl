"""
    Elastic{Dim,T<:AbstractFloat}(ρ::T, c::Complex{T})
    Elastic(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic elastic medium with wavespeed (c) and density (ρ)

Simulations in this medium produce scalar (Dim) fields in Dim dimensions. In general we use the Debye potentials to describe the field.
"""
struct Elastic{Dim,T} <: PhysicalMedium{Dim,Dim}
    ρ::T # Density (use \greekletter+tab to get the greek letter!)
    cp::Complex{T} # Phase velocity of pressure wave
    cs::Complex{T} # Phase velocity of shear wave
end

# Constructor which supplies the dimension without explicitly mentioning type
function Elastic(Dim::Integer; ρ::T = 0.0, cp::Union{T,Complex{T}} = 0.0, cs::Union{T,Complex{T}} = 0.0) where {T<:Number}
     Elastic{Dim,T}(ρ,Complex{T}(cp),Complex{T}(cs))
end

## wave types

"""
    HelmholtzPotential{Dim,T}

See [`field(::HelmholtzPotential, ::AbstractVector)`](@ref) for details on how to evaluate the Helmholtz potential.
"""
struct HelmholtzPotential{Dim,T}
    wavespeed::Complex{T}
    wavenumber::Complex{T}
    "The first (second) row is for the besselj (hankelh1) fourier coefficients"
    coefficients::Matrix{Complex{T}}
    modes::Vector{Int}

    function HelmholtzPotential{Dim}(wavespeed::Complex{T}, wavenumber::Complex{T},  coefficients::AbstractMatrix{Complex{T}}, modes::AbstractVector{Int}) where {Dim,T}
        
        # if modes |> isempty
        #     order = basislength_to_basisorder(PhysicalMedium{Dim,1},size(coefficients,2))
        #     modes = -order:order |> collect
        # end

        is = sortperm_modes(modes);

        if size(coefficients,1) != 2
            @error "the number of columns in coefficients has to match the basis_order given. There should also be two rows, one for besselj coefficients and another for hankelh1"
        end

        if imag(wavenumber) < 0
            @warn "It is usual to have a wavenumber with a negative imaginary part. Our convention of the Fourier transform implies that this wave is growing exponentially when propagating forward in time."
        end

        new{Dim,T}(wavespeed, wavenumber, coefficients[:, is], modes[is] |> collect)
    end
end


"""
    ElasticWave{Dim,T}

A type to store an elastic wave expanded in some modal basis. The only basis used in this package is a regular radial expansion, which is needed for scattering and waves within a sphere or cylinder.

For Dim = 4 the field `potentials` are, in order, the φ pressure, Φ shear,  and χ shear, where Φ and χ are the Debye potentials such that the displacement is given by: `u = ∇φ + ∇x∇x(Φ .* r) + ∇x∇x∇x(χ .* r) ./ kS`, where `kS` is the shear wavenumber and `r` is a vector pointing in the radial direction, that is `r = (x,y,z)` in Cartesian coordinates.
"""
struct ElasticWave{Dim,M,T}
    ω::T
    medium::Elastic{Dim,T}
    potentials::Vector{H} where H <:HelmholtzPotential{Dim,T}
    method::M

    function ElasticWave(ω::T, medium::Elastic{Dim,T}, potentials::Vector{H}, method::M = ModalMethod();
            mode_errors = zeros(T, length(potentials[1].modes))
        ) where {Dim,T,H <: HelmholtzPotential{Dim,T}, M<:SolutionMethod}

        modes_arr = [p.modes for p in potentials]
        lens = length.(modes_arr)

        modes = intersect(modes_arr...)
        len = length(modes)

        if !(findall(lens .!= len) |> isempty)
            @error "The modes for all the potentials is expected to be the same."
        end

        if M == ModalMethod 
            if length(method.mode_errors) != length(potentials[1].modes)
                @warn "The length of mode_errors is expected to be the same as the number of modes."
            end    
        end

        new{Dim,M,T}(ω, medium, potentials, method)
    end
end


import MultipleScattering: name, basisorder_to_basislength, basislength_to_basisorder, regular_basis_function, regular_radial_basis, outgoing_basis_function, outgoing_radial_basis, internal_field

name(a::Elastic{Dim}) where Dim = "$(Dim)D Elastic"

basisorder_to_basislength(::Type{P}, order::Int) where {T, P<:Elastic{3,T}} = 3 * (order+1)^2
basisorder_to_basislength(::Type{P}, order::Int) where {T, P<:Elastic{2,T}} = 2 * (2*order + 1)

basislength_to_basisorder(::Type{P},len::Int) where {T, P<:Elastic{3,T}} = Int(sqrt(len / 3) - 1)
basislength_to_basisorder(::Type{P},len::Int) where {T, P<:Elastic{2,T}} = Int(T(len / 2 - 1) / T(2.0))

function regular_basis_function(medium::Elastic{3,T}, ω::T, field_type::FieldType = DisplacementType()) where T

    return function (order::Integer, x::AbstractVector{T})
        pbasis = pressure_regular_basis(ω, x, medium, order, field_type) 
        Φbasis = shearΦ_regular_basis(ω, x, medium, order, field_type) 
        χbasis = shearχ_regular_basis(ω, x, medium, order, field_type) 

        # the order in each row is first iteration over p, Φ, χ, then increase basis order
        return reshape([pbasis Φbasis χbasis] |> transpose,3,:)
        
    end
end

function outgoing_basis_function(medium::Elastic{3,T}, ω::T, field_type::FieldType = DisplacementType()) where T

    return function (order::Integer, x::AbstractVector{T})
        pbasis = pressure_outgoing_basis(ω, x, medium, order, field_type)
        Φbasis = shearΦ_outgoing_basis(ω, x, medium, order, field_type)
        χbasis = shearχ_outgoing_basis(ω, x, medium, order, field_type)

        return reshape([pbasis Φbasis χbasis] |> transpose,3,:)
    end
end

"""
    internal_field(x::AbstractVector, p::Particle{Dim,Elastic{Dim,T}},  source::RegularSource{Elastic{Dim,T}}, ω::T, scattering_coefficients::AbstractVector{Complex{T}})

The internal field for an isotropic elastic particle in an isotropic elastic medium.
"""
function internal_field(x::AbstractVector{T}, p::Particle{Dim,Elastic{Dim,T}}, source::RegularSource{Elastic{Dim,T}}, ω::T, scattering_coefficients::AbstractVector{Complex{T}}, field_type::FieldType = DisplacementType()) where {Dim,T}

    if !(x ∈ p)
        @error "Point $x is not inside the particle with shape $(p.shape)"
    end
    if iszero(p.medium.cp) || isinf(abs(p.medium.cp))
        return zero(Complex{T})
    else
        fs = scattering_coefficients

        order = basislength_to_basisorder(Elastic{Dim,T},length(fs))
        r = outer_radius(p)

        t_mat = t_matrix(p, source.medium, ω, order)
        in_mat = internal_matrix(p, source.medium, ω, order)

        # need to seperate the l=0 case
        internal_coef0 = [in_mat[1] * (t_mat[1] \ fs[1]), zero(Complex{T}), zero(Complex{T})]
        
        internal_coefs = [internal_coef0; in_mat[4:end,4:end] * (t_mat[4:end,4:end] \ fs[4:end])]

        inner_basis = regular_basis_function(p.medium, ω, field_type)

        return inner_basis(order, x-origin(p)) * internal_coefs
    end
end

import Base.show
function show(io::IO, p::Elastic)
    # Print is the style of the first constructor
    write(io, "Elastic($(p.ρ), $(p.cp),  $(p.cs)) with Dim = $(spatial_dimension(p))")
    return
end

import MultipleScattering: outgoing_basis_function, regular_basis_function, outgoing_translation_matrix, regular_translation_matrix

function outgoing_basis_function(medium::Elastic{2}, ω::T) where {T<:Number}
    return function (order::Integer, x::AbstractVector{T})
        r, θ  = cartesian_to_radial_coordinates(x)
        kp = ω/medium.cp
        ks = ω/medium.cs
        vcat(
            [hankelh1(m,kp*r)*exp(im*θ*m) for m = -order:order],
            [hankelh1(m,ks*r)*exp(im*θ*m) for m = -order:order]
        ) |> transpose
    end
end

function regular_basis_function(medium::Elastic{2}, ω::T) where {T<:Number}
    return function (order::Integer, x::AbstractVector{T})
        r, θ  = cartesian_to_radial_coordinates(x)
        kp = ω/medium.cp
        ks = ω/medium.cs
        vcat(
            [besselj(m,kp*r)*exp(im*θ*m) for m = -order:order],
            [besselj(m,ks*r)*exp(im*θ*m) for m = -order:order]
        ) |> transpose
    end
end

function outgoing_translation_matrix(medium::Elastic{2}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}

    translation_vec = outgoing_basis_function(medium, ω)(in_order + out_order, x)
    order = Int(length(translation_vec)/2)
    
    translation_vec_p = translation_vec[1:order]
    translation_vec_s = translation_vec[order+1:end]
    U_p = [
        translation_vec_p[n-m + in_order + out_order + 1]
    for n in -out_order:out_order, m in -in_order:in_order]

    U_s = [
        translation_vec_s[n-m + in_order + out_order + 1]
    for n in -out_order:out_order, m in -in_order:in_order]

    U = BlockDiagonal([U_p, U_s])
    return U
end

function regular_translation_matrix(medium::Elastic{2}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}
    translation_vec = regular_basis_function(medium, ω)(in_order + out_order, x)
    order = Int(length(translation_vec)/2)
    translation_vec_p = translation_vec[1:order]
    translation_vec_s = translation_vec[order+1:end]
    V_p = [
        translation_vec_p[n-m + in_order + out_order + 1]
    for n in -out_order:out_order, m in -in_order:in_order]

    V_s = [
        translation_vec_s[n-m + in_order + out_order + 1]
    for n in -out_order:out_order, m in -in_order:in_order]

    V = BlockDiagonal([V_p, V_s])
    return V
end