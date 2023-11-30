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

import MultipleScattering: basisorder_to_basislength, basislength_to_basisorder

basisorder_to_basislength(::Type{P}, order::Int) where {T, P<:Elastic{T,3}} = 3 * (order+1)^2
basisorder_to_basislength(::Type{P}, order::Int) where {T, P<:Elastic{T,2}} = 2 * (2*order + 1)

basislength_to_basisorder(::Type{P},len::Int) where {T, P<:Elastic{T,3}} = Int(sqrt(len / 3) - 1)
basislength_to_basisorder(::Type{P},len::Int) where {T, P<:Elastic{T,2}} = Int(T(len / 2 - 1) / T(2.0))

import Base.show
function show(io::IO, p::Elastic)
    # Print is the style of the first constructor
    write(io, "Elastic($(p.ρ), $(p.cp),  $(p.cs)) with Dim = $(spatial_dimension(p))")
    return
end

name(a::Elastic{Dim}) where Dim = "$(Dim)D Elastic"


## field types
"""
    FieldType

A type used to specify what type of physical field, such as traction or displacement.
"""
abstract type FieldType end

struct DisplacementType <: FieldType end
struct TractionType <: FieldType end

## wave types

"""
    HelmholtzPotential{Dim,T}

See [`field(::HelmholtzPotential, ::AbstractVector)`](@ref) for details on how to evaluate the Helmholtz potential.
"""
struct HelmholtzPotential{Dim,T}
    wavespeed::Complex{T}
    wavenumber::Complex{T}
    basis_order::Int
    "The first (second) row is for the besselj (hankelh1) fourier modes"
    coefficients::Matrix{Complex{T}}

    function HelmholtzPotential{Dim}(wavespeed::Complex{T}, wavenumber::Complex{T},  coefficients::AbstractMatrix{Complex{T}},
            basis_order::Int = basislength_to_basisorder(Acoustic{T,Dim},size(coefficients,2))
        ) where {Dim,T}
        if size(coefficients,1) != 2
            @error "the number of columns in coefficients has to match the basis_order given. There should also be two rows, one for besselj coefficients and another for hankelh1"
        end

        if imag(wavenumber) < 0
            @warn "It is usual to have a wavenumber with a negative imaginary part. Our convention of the Fourier transform implies that this wave is growing exponentially when propagating forward in time."
        end

        new{Dim,T}(wavespeed, wavenumber, basis_order, coefficients)
    end
end

struct ElasticWave{Dim,T}
    ω::T
    medium::Elastic{Dim,T}
    pressure::HelmholtzPotential{Dim,T}
    shear::HelmholtzPotential{Dim,T}
    mode_errors::Vector{T}

    function ElasticWave(ω::T, medium::Elastic{Dim,T}, pressure::HelmholtzPotential{Dim,T}, shear::HelmholtzPotential{Dim,T};
            mode_errors = zeros(basisorder_to_basislength(Acoustic{T,Dim}, pressure.basis_order))
        ) where {Dim,T}

        if pressure.basis_order != shear.basis_order
            @error "The basis_order for both potentials is expected to be the same."
        end

        if length(mode_errors) != basisorder_to_basislength(Acoustic{T,Dim}, pressure.basis_order)
            @error "The length of mode_errors should be the same as the number of modes."
        end

        new{Dim,T}(ω, medium, pressure, shear, mode_errors)
    end
end
