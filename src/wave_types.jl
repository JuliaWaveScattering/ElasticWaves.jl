"""
    Elasticity{Dim,T<:AbstractFloat}(ρ::T, c::Complex{T})
    Elasticity(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic elastic medium with wavespeed (c) and density (ρ)

Simulations in this medium produce scalar (Dim) fields in Dim dimensions.
"""
struct Elasticity{Dim,T} <: PhysicalMedium{Dim,Dim}
    ρ::T # Density (use \greekletter+tab to get the greek letter!)
    cp::Complex{T} # Phase velocity of pressure wave
    cs::Complex{T} # Phase velocity of shear wave
end


# Constructor which supplies the dimension without explicitly mentioning type
function Elasticity(Dim::Integer; ρ::T = 0.0, cp::Union{T,Complex{T}} = 0.0, cs::Union{T,Complex{T}} = 0.0) where {T<:Number}
     Elasticity{Dim,T}(ρ,Complex{T}(cp),Complex{T}(cs))
end

import Base.show
function show(io::IO, p::Elasticity)
    # Print is the style of the first constructor
    write(io, "Elasticity($(p.ρ), $(p.cp),  $(p.cs)) with Dim = $(spatial_dimension(p))")
    return
end

name(a::Elasticity{Dim}) where Dim = "$(Dim)D Elasticity"

## wave types

"""
    HelmholtzPotential{Dim,T}

See [`field(::HelmholtzPotential, ::AbstractVector)`](@ref) for details on how to evaluate the Helmholtz potential.
"""
struct HelmholtzPotential{Dim,T}
    wavenumber::Complex{T}
    basis_order::Int
    "The first (second) row is for the besselj (hankelh1) fourier modes"
    coefficients::Matrix{Complex{T}}

    function HelmholtzPotential{Dim}(wavenumber::Complex{T}, coefficients::AbstractMatrix{Complex{T}},
            basis_order::Int = basislength_to_basisorder(Acoustic{T,Dim},size(coefficients,2))
        ) where {Dim,T}
        if size(coefficients,1) != 2
            @error "the number of rows in coefficients has to match the basis_order given. There should also be two columns, one for besselj coefficients and another for hankelh1"
        end

        if imag(wavenumber) < 0
            @warn "It is usual to have a wavenumber with a negative imaginary part. Our convention of the Fourier transform implies that this wave is growing exponentially when propagating forward in time."
        end

        new{Dim,T}(wavenumber, basis_order, coefficients)
    end
end

struct ElasticWave{Dim,T}
    ω::T
    medium::Elasticity{Dim,T}
    pressure::HelmholtzPotential{Dim,T}
    shear::HelmholtzPotential{Dim,T}

    function ElasticWave{Dim}(ω::T, medium::Elasticity{Dim,T}, pressure::HelmholtzPotential{Dim,T}, shear::HelmholtzPotential{Dim,T}) where {Dim,T}
        if pressure.basis_order != shear.basis_order
            @error "The basis_order for both potentials is expected to be the same."
        end
        new{Dim,T}(ω, medium,pressure,shear)
    end
end
