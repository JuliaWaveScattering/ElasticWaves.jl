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

    println("Note that the traction on the inner boundary (τn,τt) is interpreted to be the vector τ = - τn * er - τt * eθ.")

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

    function BoundaryBasis(basis::Vector{BoundaryData{BC,T}}) where {BC <: BoundaryCondition, T}
        field_lengths = [size(b.fields, 1) for b in basis]
        fourier_mode_lengths = [size(b.fourier_modes, 1) for b in basis]

        allsame(x) = all(y -> y == first(x), x)

        if !allsame(fourier_mode_lengths) || !allsame(field_lengths)
            error("The length of the fields and fourier_modes of each BoundaryData in a BoundaryBasis need to be equal")
        end    

    return new{BC,T}(basis)
    end    
end

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

mutable struct BearingSimulation{M <: SolutionMethod, BC1 <: BoundaryCondition, BC2 <: BoundaryCondition, BCB1 <: BoundaryCondition, BCB2 <: BoundaryCondition, T}
    ω::T
    basis_order::Int
    method::M
    bearing::RollerBearing{T}
    boundarydata1::BoundaryData{BC1,T}
    boundarydata2::BoundaryData{BC2,T}
    boundarybasis1::BoundaryBasis{BCB1,T}
    boundarybasis2::BoundaryBasis{BCB2,T}  
    
    function BearingSimulation(ω::T, basis_order::Int, method::M, bearing::RollerBearing{T}, 
            boundarydata1::BoundaryData{BC1,T},
            boundarydata2::BoundaryData{BC2,T},
            boundarybasis1::BoundaryBasis{BCB1,T}  = BoundaryBasis([BoundaryData(boundarydata1.boundarytype)]),
            boundarybasis2::BoundaryBasis{BCB2,T}  = BoundaryBasis([BoundaryData(boundarydata2.boundarytype)])
        ) where {T, M <: SolutionMethod, BC1 <: BoundaryCondition, BC2 <: BoundaryCondition, BCB1 <: BoundaryCondition, BCB2 <: BoundaryCondition}
    
        return new{M,BC1,BC2,BCB1,BCB2,T}(ω, basis_order, method, bearing, boundarydata1, boundarydata2, boundarybasis1, boundarybasis2)
    end 
end


function BearingSimulation(ω::T, bearing::RollerBearing{T}, boundarydata1::BD1, boundarydata2::BD2;
    method::M = NoBearingMethod(),
    kws...
) where {T, M <: SolutionMethod, BD1 <: BoundaryData, BD2 <: BoundaryData}

    if typeof(method) == NoBearingMethod
        @warn "no method was chosen for this bearing simulation. Will choose the default ModalMethod which makes no assumption about the forcing or bearings"
        
        method = ModalMethod()
    end    
    
    return BearingSimulation(ω ,method, bearing, boundarydata1, boundarydata2; kws...)
end

function BearingSimulation(ω::T, method::ModalMethod, bearing::RollerBearing{T}, boundarydata1::BoundaryData{BC1,T}, boundarydata2::BoundaryData{BC2,T};
        basis_order::Int = -1,
    ) where {T, BC1 <: BoundaryCondition, BC2 <: BoundaryCondition}

    if basis_order == -1
        if isempty(boundarydata1.fourier_modes)
            if isempty(boundarydata2.fourier_modes)
                println("As the keyword basis_order was not specified (and the fourier_modes of the boundary conditions were not provided), the basis_order will be chosen to match the number of points θs given in the boundary data.")

                # lower the basis_order if there are not enough data points to estimate it
                m1 = Int(floor(length(boundarydata1.θs)/2.0 - 1/2.0))
                m2 = Int(floor(length(boundarydata2.θs)/2.0 - 1/2.0))

                basis_order = min(m1, m2)

            else basis_order = basislength_to_basisorder(Acoustic{T,2},size(boundarydata2.fourier_modes,1))
            end

        else basis_order = basislength_to_basisorder(Acoustic{T,2},size(boundarydata1.fourier_modes,1))
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

    if size(boundarydata1.fourier_modes,1) != size(boundarydata2.fourier_modes,1)
        error("number of fourier_modes in boundarydata1 and boundarydata2 needs to be the same")
    end

    return BearingSimulation(ω, basis_order, method, bearing, boundarydata1, boundarydata2)
end

function BearingSimulation(ω::T, method::PriorMethod, bearing::RollerBearing{T}, boundarydata1::BoundaryData{BC1,T}, boundarydata2::BoundaryData{BC2,T};
        basis_order::Int = -1,
        boundarybasis1::BoundaryBasis{BCB1} = BoundaryBasis([BoundaryData(boundarydata1.boundarytype)]),
        boundarybasis2::BoundaryBasis{BCB2} = BoundaryBasis([BoundaryData(boundarydata2.boundarytype)])
    ) where {T, BC1 <: BoundaryCondition, BC2 <: BoundaryCondition, BCB1 <: BoundaryCondition, BCB2 <: BoundaryCondition}

    # we need the fields of the boundary data for the PriorMethod
    if isempty(boundarydata1.fields) && isempty(boundarydata2.fields)

        @error "The fields of boundarydata1 or boundarydata2 are empty. We expect that the fields are given for the PriorMethod because as the number of fourier_modes can be far greater than the length of the fields. Meaning it is not possible to calculate the fourier_modes of the boundary data directly"
    end

    if basis_order == -1
        println("No basis_order was specified. Will use the largest order according to the data provided for boundarybasis1 and boundarybasis2")

        # if both boundarybasis1 and boundarybasis2 are given need to use the same number of modes
        mode_lengths = [size(bb.basis[1].fourier_modes,1) for bb in [boundarybasis1, boundarybasis2]]

        inds = findall(mode_lengths .> 0)
        mode_lengths = unique(mode_lengths[inds])

        if length(mode_lengths) == 2
            error("The number of Fourier modes for boundarybasis1 and boundarybasis2 need to be the same")
        end
        
        if isempty(mode_lengths)
            field_lengths = [size(bb.basis[1].fields,1) for bb in [boundarybasis1, boundarybasis2]]
            inds = findall(field_lengths .> 0)
            field_lengths = field_lengths[inds]

            # need to use the smallest field length to calculate the fourier_modes of both boundarybasis'
            field_length = minimum(field_lengths)
            basis_order = Int(floor((field_length - 1)/2 ))

        else basis_order = basislength_to_basisorder(PhysicalMedium{2,1},mode_lengths[1])
        end
    end

    # We need the fourier_modes for the boundarybasis
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

    return BearingSimulation(ω, basis_order, method, bearing, boundarydata1, boundarydata2, boundarybasis1, boundarybasis2)
end
