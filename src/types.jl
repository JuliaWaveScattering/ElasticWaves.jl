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

struct HelmholtzPotential{Dim,T}
    wavenumber::Complex{T}
    basis_order::Int
    "The first (second) row is for besselj (hankelh1) modes"
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

    function ElasticWave{Dim}(ω::T, medium::Elasticity{Dim,T},pressure::HelmholtzPotential{Dim,T}, shear::HelmholtzPotential{Dim,T}) where {Dim,T}
        if pressure.basis_order != shear.basis_order
            @error "The basis_order for both potentials is expected to be the same."
        end
        new{Dim,T}(ω, medium,pressure,shear)
    end
end


struct RollerBearing{T}
    medium::Elasticity{2,T} # defining medium
    r1::T # inner radius
    r2::T # outer radius
end

function RollerBearing(; medium::Elasticity{2}, r1::T=0.0, r2::T=0.0) where {T<:Number}
    RollerBearing{T}(medium,r1,r2)
end


## Boundary condtion types

abstract type AbstractBoundaryCondition end

AbstractBoundaryConditions = Vector{BC} where BC <: AbstractBoundaryCondition

# Art NOTE: we should create a type called BearingSystem, which holds both the forcing and the BCs.

struct DisplacementBoundary <: AbstractBoundaryCondition
    inner::Bool
    outer::Bool
end

function DisplacementBoundary(; inner::Bool = false, outer::Bool = false) where T
    if inner == outer
        @error "this type represents only one boundary, either the inner or outer boundary, and not both or neither."
    end

    return DisplacementBoundary(inner,outer)
end

struct TractionBoundary <: AbstractBoundaryCondition
    inner::Bool
    outer::Bool
end

function TractionBoundary(; inner::Bool = false, outer::Bool = false) where T
    if inner == outer
        @error "this type represents only one boundary, either the inner or outer boundary, and not both or neither."
    end

    return TractionBoundary(inner,outer)
end


"""
    BoundaryData{BC <: AbstractBoundaryCondition, T}

Gives all the information on one boundery of a `RollerBearing`. You can specify the forcing either be giving the `modes`, or by specifying the `fields`. The `fields` represent either the displacement or traction when varying the polar angle `θs`.

The relationship between the `fields` and the modes is given by ``f_{j}(\\theta) = \\sum_m \\text{mode}_{mj} e^{i \\theta m}``, where ``f`` represents the `fields`. To calculate the `fields` from the `modes` we can use:
```julia
basis_order = (size(modes,1) - 1) / 2
exps = [
    exp(im * θ * m)
for θ in θs, m = -basis_order:basis_order];

field = exps * modes
```

If the `modes` were specified, then they will be used directly. If the `modes` were not give, then `fields` will be used to calculate the `modes`: `modes = exps \\ fields`.
"""
struct BoundaryData{BC <: AbstractBoundaryCondition, T}
    boundarytype::BC
    θs::Vector{T}
    fields::Matrix{Complex{T}}
    modes::Matrix{Complex{T}}
end

function BoundaryData(boundarytype::BC, θs::Vector{T} = Float64[], fields::Matrix = reshape(Complex{Float64}[],0,1);
        modes::Matrix = reshape(Complex{Float64}[],0,1)
    ) where {BC <: AbstractBoundaryCondition, T}

    if !isempty(modes) && size(modes,2) != 2
        @error "the modes should have only two columns"
    end

    if !isempty(fields) && size(fields,2) != 2
        @error "the fields should have only two columns"
    end

    return BoundaryData{BC,T}(boundarytype,θs,fields,modes)
end


struct BearingSimulation{BC1 <: AbstractBoundaryCondition, BC2 <: AbstractBoundaryCondition, T}
    ω::T
    basis_order::Int
    bearing::RollerBearing{T}
    boundarydata1::BoundaryData{BC1,T}
    boundarydata2::BoundaryData{BC2,T}
    tol::T

    function BearingSimulation(ω::T, bearing::RollerBearing{T}, boundarydata1::BoundaryData{BC1,T}, boundarydata2::BoundaryData{BC2,T}; tol::T = eps(T)^(1/2)) where {T, BC1 <: AbstractBoundaryCondition, BC2 <: AbstractBoundaryCondition}

        if size(boundarydata1.modes,1) != size(boundarydata2.modes,1)
            @error "number of modes in boundarydata1 and boundarydata2 needs to be the same"
        end

        basisorder = basislength_to_basisorder(Acoustic{T,2},size(boundarydata1.modes,1))

        new{BC1,BC2,T}(ω, basisorder, bearing, boundarydata1, boundarydata2, tol)
    end
end
