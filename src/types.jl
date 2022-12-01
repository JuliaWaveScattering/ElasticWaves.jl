"""
    Elasticity{T<:AbstractFloat,Dim}(ρ::T, c::Complex{T})
    Elasticity(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic elastic medium with wavespeed (c) and density (ρ)

Simulations in this medium produce scalar (Dim) fields in Dim dimensions.
"""
struct Elasticity{T,Dim} <: PhysicalMedium{Dim,Dim}
    ρ::T # Density (use \greekletter+tab to get the greek letter!)
    cp::Complex{T} # Phase velocity of pressure wave
    cs::Complex{T} # Phase velocity of shear wave
end


# Constructor which supplies the dimension without explicitly mentioning type
function Elasticity(Dim::Integer; ρ::T = 0.0, cp::Union{T,Complex{T}} = 0.0, cs::Union{T,Complex{T}} = 0.0) where {T<:Number}
     Elasticity{T,Dim}(ρ,Complex{T}(cp),Complex{T}(cs))
end

import Base.show
function show(io::IO, p::Elasticity)
    # Print is the style of the first constructor
    write(io, "Elasticity($(p.ρ), $(p.cp),  $(p.cs)) with Dim = $(spatial_dimension(p))")
    return
end

name(a::Elasticity{T,Dim}) where {Dim,T} = "$(Dim)D Elasticity"

struct Bearing{T,Dim}
    medium::Elasticity{T,Dim} # defining medium
    r1::T # inner radius
    r2::T # outer radius
end

function Bearing(Dim::Int; medium::Elasticity, r1::T=0.0,r2::T=0.0) where {T<:Number}
    Bearing{T,Dim}(medium,r1,r2)
end

## Boundary condtion types

abstract type AbstractBoundaryCondition end

AbstractBoundaryConditions = Vector{BC} where BC <: AbstractBoundaryCondition

# Art NOTE: we should create a type called BearingBoundaryConditions, which holds both the forcing and the BCs.

struct DisplacementBoundary <: AbstractBoundaryCondition
    inner::Bool
    outer::Bool
end

function DisplacementBoundary(;inner::Bool = false, outer::Bool = false)
    if inner == outer
        @error "this type represents only one boundary, either the inner or outer boundary, and not both or neither."
    end

    return DisplacementBoundary(inner,outer)
end

struct StressBoundary <: AbstractBoundaryCondition
    inner::Bool
    outer::Bool
end

function StressBoundary(;inner::Bool = false, outer::Bool = false)
    if inner == outer
        @error "this type represents only one boundary, either the inner or outer boundary, and not both or neither."
    end

    return StressBoundary(inner,outer)
end
