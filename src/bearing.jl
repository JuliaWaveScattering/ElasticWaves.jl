## methods to solve for waves in bearings 

abstract type BearingMethod end

struct ModalMethod <: BearingMethod end
struct GapMethod <: BearingMethod end
struct PriorMethod <: BearingMethod end


struct RollerBearing{T}
    medium::Elastic{2,T} # defining medium
    inner_radius::T
    "vector of angles delimiting gaps in the inner radius"
    inner_gaps::Vector{T}
    outer_radius::T
    "vector of angles delimiting gaps in the outer radius"
    outer_gaps::Vector{T}
    number_of_rollers::T
end

function RollerBearing(; medium::Elastic{2},
        inner_radius::T = 0.0, outer_radius::T = 0.0,
        inner_gaps::Vector{T} = typeof(inner_radius)[],
        outer_gaps::Vector{T} = typeof(outer_radius)[],
        number_of_rollers::T=1.0
    ) where T<:Number
    if isodd(length(inner_gaps)) && isodd(length(outer_gaps))
        @error "both inner_gaps and outer_gaps need to be an even number of angles"
    end
    RollerBearing{T}(medium, inner_radius, inner_gaps, outer_radius, outer_gaps, number_of_rollers)
end


## a type to represent one of the Boundary condtions

struct BoundaryCondition{FT<:FieldType}
    fieldtype::FT
    inner::Bool
    outer::Bool
end

BoundaryConditions = Vector{BC} where BC <: BoundaryCondition

# Art NOTE: we should create a type called BearingSystem, which holds both the forcing and the BCs.

function DisplacementBoundary(; inner::Bool = false, outer::Bool = false)
    if inner == outer
        @error "this type represents only one boundary, either the inner or outer boundary, and not both or neither."
    end

    return BoundaryCondition{DisplacementType}(DisplacementType(),inner,outer)
end

function TractionBoundary(; inner::Bool = false, outer::Bool = false)
    if inner == outer
        @error "this type represents only one boundary, either the inner or outer boundary, and not both or neither."
    end

    return BoundaryCondition{TractionType}(TractionType(),inner,outer)
end


"""
    BoundaryData{BC <: BoundaryCondition, T}

Gives all the information on one boundery of a `RollerBearing`. You can specify the forcing either be giving the `fourier_modes`, or by specifying the `fields`. The `fields` represent either the displacement or traction when varying the polar angle `θs`.

The relationship between the `fields` and the `fourier_modes` is given by ``f_{j}(\\theta) = \\sum_m \\text{mode}_{mj} e^{i \\theta m}``, where ``f`` represents the `fields`. To calculate the `fields` from the `fourier_modes` we can use:
```julia
basis_order = (size(fourier_modes,1) - 1) / 2
exps = [
    exp(im * θ * m)
for θ in θs, m = -basis_order:basis_order];

fields = exps * fourier_modes
```

If the `fourier_modes` were specified, then they will be used directly. If the `fourier_modes` were not give, then `fields` will be used to calculate the `fourier_modes`: `fourier_modes = exps \\ fields`.
"""
struct BoundaryData{BC <: BoundaryCondition, T}
    boundarytype::BC
    θs::Vector{T}
    fields::Matrix{Complex{T}}
    fourier_modes::Matrix{Complex{T}}
end

function BoundaryData(boundarytype::BC;
        θs::AbstractVector{T} = Float64[],
        fields::Matrix = reshape(Complex{Float64}[],0,2),
        fourier_modes::Matrix = reshape(Complex{Float64}[],0,2)
    ) where {BC <: BoundaryCondition, T}

    if !isempty(fourier_modes) && size(fourier_modes,2) != 2
        @error "the fourier_modes should have only two columns"
    end

    if !isempty(fields) && size(fields,2) != 2
        @error "the fields should have only two columns"
    end

    return BoundaryData{BC,T}(boundarytype,θs,fields,fourier_modes)
end

struct BoundaryBasis{BC <: BoundaryCondition, T}
    basis::Vector{BoundaryData{BC,T}}
end

# function BoundaryBasis(basis::Vector{BoundaryData{BC,T}}) where {BC <: BoundaryCondition, T}
#     return BoundaryBasis{BC, T}(basis)
# end

# """
#     EmptyBoundaryData(T)

# A shorthand for what will be considered an empty  BoundaryData
# """
# EmptyBoundaryData() =  BoundaryData(BoundaryCondition{DisplacementType}(DisplacementType(),false,false))

# """
#     EmptyBoundaryBasis(T)

# A shorthand for what will be considered an empty BoundaryBasis
# """
# EmptyBoundaryBasis(bd::BoundaryData) =  BoundaryBasis(bd[])

struct BearingSimulation{M <: BearingMethod, BC1 <: BoundaryCondition, BC2 <: BoundaryCondition, BCB1 <: BoundaryCondition, BCB2 <: BoundaryCondition, T}
    ω::T
    basis_order::Int
    method::M
    tol::T
    bearing::RollerBearing{T}
    boundarydata1::BoundaryData{BC1,T}
    boundarydata2::BoundaryData{BC2,T}
    boundarybasis1::BoundaryBasis{BCB1,T}
    boundarybasis2::BoundaryBasis{BCB2,T}
end

function BearingSimulation(ω::T, bearing::RollerBearing{T}, boundarydata1::BoundaryData{BC1,T}, boundarydata2::BoundaryData{BC2,T};
        tol::T = eps(T)^(1/2), basis_order::Int = -1,
        boundarybasis1::BoundaryBasis{BCB1} = BoundaryBasis([BoundaryData(boundarydata1.boundarytype)]),
        boundarybasis2::BoundaryBasis{BCB2} = BoundaryBasis([BoundaryData(boundarydata2.boundarytype)])
    ) where {T, M <: BearingMethod, BC1 <: BoundaryCondition, BC2 <: BoundaryCondition, BCB1 <: BoundaryCondition, BCB2 <: BoundaryCondition}

    if size(boundarydata1.fourier_modes,1) != size(boundarydata2.fourier_modes,1)
        @error "number of fourier_modes in boundarydata1 and boundarydata2 needs to be the same"
    end

    if basis_order == -1
        if isempty(boundarydata1.fourier_modes)
            if isempty(boundarydata2.fourier_modes)
                println("As the keyword basis_order was not specified (and the fourier_modes of the boundary conditions were not provided), the basis_order will be estimated from the bearing geometry and wavenumbers.")

                basis_order = estimate_basisorder(ω, bearing; tol = 1e3 * tol)
                # lower the basis_order if there are not enough data points to estimate it
                m1 = Int(round(length(boundarydata1.θs)/2.0 - 1/2.0))
                m2 = Int(round(length(boundarydata2.θs)/2.0 - 1/2.0))

                basis_order = min(m1, m2, basis_order)

            else basis_order = basislength_to_basisorder(Acoustic{T,2},size(boundarydata2.fourier_modes,1))
            end

        else basis_order = basislength_to_basisorder(Acoustic{T,2},size(boundarydata2.fourier_modes,1))
        end
    end

   # we need the fourier modes of the boundary data
    if isempty(boundarydata1.fourier_modes)
        println("The Fourier modes for boundarydata1 are empty, they will be calculated from the fields provided")

        boundarydata1 = fields_to_fouriermodes(boundarydata1, basis_order)
    end

    if isempty(boundarydata2.fourier_modes)
        println("The Fourier modes for boundarydata2 are empty, they will be calculated from the fields provided")

        boundarydata2 = fields_to_fouriermodes(boundarydata2, basis_order)
    end

    # Calculate the fourier_modes for the boundarybasis
    basis_vec = map(boundarybasis1.basis) do b
        if isempty(b.fourier_modes)
            isempty(b.fields) ? b : fields_to_fouriermodes(b,basis_order)
        else
            b
        end
    end
    boundarybasis1 = BoundaryBasis(basis_vec)
    
    basis_vec = map(boundarybasis2.basis) do b
        if isempty(b.fourier_modes)
            isempty(b.fields) ? b : fields_to_fouriermodes(b,basis_order)
        else
            b
        end
    end
    boundarybasis2 = BoundaryBasis(basis_vec)

    if !isempty(bearing.inner_gaps) || !isempty(bearing.outer_gaps)
        method = GapMethod()
    elseif !isempty(boundarybasis1.basis[1].fourier_modes) || !isempty(boundarybasis2.basis[1].fourier_modes)
        method = PriorMethod()
    else 
        method = ModalMethod()   
    end

    
    
    BearingSimulation{typeof(method),BC1,BC2,BCB1,BCB2,T}(ω, basis_order, method, tol, bearing, boundarydata1, boundarydata2, boundarybasis1, boundarybasis2)
end
