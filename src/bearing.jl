struct RollerBearing{T <: AbstractFloat}
    medium::Elastic{2,T} # defining medium
    inner_radius::Union{T,Complex{T}} # due to nondimensionalisation radius might be complex
    "vector of angles delimiting gaps in the inner radius"
    inner_gaps::Vector{T}
    outer_radius::Union{T,Complex{T}}
    "vector of angles delimiting gaps in the outer radius"
    outer_gaps::Vector{T}
    number_of_rollers::Int
    rollers_inside::Boole
    angular_speed::Float64
end

function RollerBearing(; medium::Elastic{2,T},
        inner_radius::Union{T,Complex{T}} = 0.0, 
        outer_radius::Union{T,Complex{T}} = 0.0,
        inner_gaps::Vector{T} = typeof(medium).parameters[2][],
        outer_gaps::Vector{T} = typeof(medium).parameters[2][],
        number_of_rollers::Int = 1,
        rollers_inside = true,
        angular_speed=1.0
    ) where T <: AbstractFloat
    if isodd(length(inner_gaps)) && isodd(length(outer_gaps))
        @error "both inner_gaps and outer_gaps need to be an even number of angles"
    end

    println("Note that the traction on the inner boundary (τn,τt) is interpreted to be the vector τ = - τn * er - τt * eθ.")

    RollerBearing{T}(medium, inner_radius, inner_gaps, outer_radius, outer_gaps, number_of_rollers,rollers_inside,angular_speed)
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

abstract type AbstractBoundaryData{BC <: BoundaryCondition} end   

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
struct BoundaryData{BC <: BoundaryCondition, T} <: AbstractBoundaryData{BC}
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

struct LoadingProfile{BC <: BoundaryCondition, T} <: AbstractBoundaryData{BC}
    boundarytype::BC
    θs::Vector{T}
    fields::Matrix{Complex{T}}
    fourier_modes::Matrix{Complex{T}}
end

function LoadingProfile(boundarytype::BC;
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

import Base: isempty
isempty(bd::BoundaryData) = all([isempty(bd.fields),isempty(bd.fourier_modes)])  
isempty(bb::BoundaryBasis) = all([isempty(bd) for bd in bb.basis])  
    
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
    method::M
    nondimensionalise::Bool 
    bearing::RollerBearing{T}
    boundarydata1::BoundaryData{BC1,T}
    boundarydata2::BoundaryData{BC2,T}
    boundarybasis1::BoundaryBasis{BCB1,T}
    boundarybasis2::BoundaryBasis{BCB2,T}  
end

function BearingSimulation(ω::T, method::M, bearing::RollerBearing{T}, 
        boundarydata1::BoundaryData{BC1,T},
        boundarydata2::BoundaryData{BC2,T};
        nondimensionalise::Bool = true,
        boundarybasis1::BoundaryBasis{BCB1,T}  = BoundaryBasis([BoundaryData(boundarydata1.boundarytype)]),
        boundarybasis2::BoundaryBasis{BCB2,T}  = BoundaryBasis([BoundaryData(boundarydata2.boundarytype)])
    ) where {T, M <: SolutionMethod, BC1 <: BoundaryCondition, BC2 <: BoundaryCondition, BCB1 <: BoundaryCondition, BCB2 <: BoundaryCondition}

    sim = BearingSimulation{M,BC1,BC2,BCB1,BCB2,T}(ω, method, nondimensionalise, bearing, boundarydata1, boundarydata2, boundarybasis1, boundarybasis2)

    return setup(sim)
end 

import MultipleScattering: boundary_data
"""
    boundary_data(wave::ElasticWave{2}, sim::BearingSimulation)

A convenience function to predict the boundary data of 'sim' according to the solution 'wave'.    
"""
function boundary_data(sim::BearingSimulation, wave::ElasticWave{2})

    bd1 = sim.boundarydata1;
    bd2 = sim.boundarydata2;

    radius_inner_outer(in) = (in == true) ? sim.bearing.inner_radius : sim.bearing.outer_radius

    inners = [bd1.boundarytype.inner, bd2.boundarytype.inner]
    radiuses = radius_inner_outer.(inners)
    
    bd1 = BoundaryData(bd1.boundarytype, radiuses[1], bd1.θs, wave)
    bd2 = BoundaryData(bd2.boundarytype, radiuses[2], bd2.θs, wave)

    return [bd1, bd2]
end

"""
    BoundaryData(wave::ElasticWave{2}, radius, θs, fieldtype)

A convenience function to predict the boundary data of 'sim' according to the solution 'wave'.    
"""
function BoundaryData(boundarycondition::BoundaryCondition, radius::AbstractFloat, θs::AbstractVector{T}, wave::ElasticWave{2}) where T <: AbstractFloat

    xs = [
        radial_to_cartesian_coordinates([radius,θ]) 
    for θ in θs];
    
    fs = [field(wave, x, boundarycondition.fieldtype) for x in xs];
    fs = hcat(fs...) |> transpose |> collect

    return BoundaryData(boundarycondition; θs = θs, fields = fs)
end

function nondimensionalise!(boundarydata::BoundaryData{BoundaryCondition{TractionType}}, sim::BearingSimulation)

    λ2μ = sim.bearing.medium.ρ * sim.bearing.medium.cp^2
    boundarydata.fields[:] = boundarydata.fields ./ λ2μ
    boundarydata.fourier_modes[:] = boundarydata.fourier_modes ./ λ2μ

    return boundarydata
end

function nondimensionalise!(boundarydata::BoundaryData{BoundaryCondition{DisplacementType}}, sim::BearingSimulation)

    kp = sim.ω / sim.bearing.medium.cp;
    boundarydata.fields[:] = boundarydata.fields .* kp
    boundarydata.fourier_modes[:] = boundarydata.fourier_modes .* kp

    return boundarydata
end

function nondimensionalise(bearing::RollerBearing, sim::BearingSimulation)

    kp = sim.ω / bearing.medium.cp;
    non_medium = Elastic(2; ρ = 1 / sim.ω^2, cp = sim.ω, cs = bearing.medium.cs * kp)
    
    # nondimensionalise bearing geometry
    @reset bearing.medium = non_medium
    @reset bearing.inner_radius = bearing.inner_radius * kp
    @reset bearing.outer_radius = bearing.outer_radius * kp
    
    return bearing
end

function nondimensionalise!(sim::BearingSimulation)

    # nondimensionalise boundary data
    nondimensionalise!(sim.boundarydata1, sim)
    nondimensionalise!(sim.boundarydata2, sim)
    
    # nondimensionalise boundary basis
    sim.boundarybasis1.basis[:] = [nondimensionalise!(bd, sim) for bd in sim.boundarybasis1.basis]
    sim.boundarybasis2.basis[:] = [nondimensionalise!(bd, sim) for bd in sim.boundarybasis2.basis]

    # nondimensionalise bearing
    sim.bearing = nondimensionalise(sim.bearing,sim)

    return sim
end


function BearingSimulation(ω::T, bearing::RollerBearing{T}, boundarydata1::BD1, boundarydata2::BD2;
    method::M = NoBearingMethod(),
    kws...
) where {T, M <: SolutionMethod, BD1 <: BoundaryData, BD2 <: BoundaryData}

    if typeof(method) == NoBearingMethod
        @warn "no method was chosen for this bearing simulation. Will choose the default ModalMethod which makes no assumption about the forcing or bearings"
        
        method = ModalMethod()
    end    
    
    return BearingSimulation(ω, method, bearing, boundarydata1, boundarydata2; kws...) 
end

function setup(sim::BearingSimulation{ModalMethod})
    
    basis_order = sim.method.basis_order
    boundarydata1 = sim.boundarydata1
    boundarydata2 = sim.boundarydata2

    if basis_order == -1
        if isempty(boundarydata1.fourier_modes)
            if isempty(boundarydata2.fourier_modes)
                println("As the keyword basis_order was not specified (and the fourier_modes of the boundary conditions were not provided), the basis_order will be chosen to match the number of points θs given in the boundary data.")

                # lower the basis_order if there are not enough data points to estimate it
                m1 = Int(floor(length(boundarydata1.θs)/2.0 - 1/2.0))
                m2 = Int(floor(length(boundarydata2.θs)/2.0 - 1/2.0))

                basis_order = min(m1, m2)

            else basis_order = basislength_to_basisorder(PhysicalMedium{2,1},size(boundarydata2.fourier_modes,1))
            end

        else basis_order = basislength_to_basisorder(PhysicalMedium{2,1},size(boundarydata1.fourier_modes,1))
        end

        @reset sim.method.basis_order = basis_order
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

    @reset sim.boundarydata1 = boundarydata1
    @reset sim.boundarydata2 = boundarydata2

    return sim
end

function setup(sim::BearingSimulation{PriorMethod})

    basis_order = sim.method.modal_method.basis_order;

    if isempty(sim.boundarybasis1) && isempty(sim.boundarybasis2)
        error("The PriorMethod requires a boundary basis either for the first boundary (keyword 'boundary_basis1') or the second boundary (keyword 'boundary_basis2') ")
    end    

    # we need the fields of the boundary data for the PriorMethod
    if isempty(sim.boundarydata1.fields) && isempty(sim.boundarydata2.fields)

        error("The fields of boundarydata1 or boundarydata2 are empty. We expect that the fields are given for the PriorMethod because as the number of fourier_modes can be far greater than the length of the fields. Meaning it is not possible to calculate the fourier_modes of the boundary data directly")
    end

    if basis_order == -1
        println("No basis_order was specified. Will use the largest order according to the data provided (boundarybasis1, boundarybasis2, and poentially the boundarydata1 or boundarydata2)")

        # need the Fourier modes for two boundaries. This can be provided either by one BoundaryBasis and one BoundaryData (most common), or by two BoundaryBasis. We determine the scenario first:

        data1 = if isempty(sim.boundarybasis1)
            sim.boundarydata1
        else
            sim.boundarybasis1.basis[1]
        end
        
        data2 = if isempty(sim.boundarybasis2)
            sim.boundarydata2
        else
            sim.boundarybasis2.basis[1]
        end

        # if both boundarybasis1 and boundarybasis2 are given need to use the same number of modes
        mode_lengths = [size(bd.fourier_modes,1) for bd in [data1, data2]]

        inds = findall(mode_lengths .> 0)
        mode_lengths = unique(mode_lengths[inds])

        if length(mode_lengths) == 2
            warn("The number of Fourier modes for the data need to be the same but were different.")
            mode_lengths = minimum(mode_lengths)
        end
        
        if isempty(mode_lengths)
            field_lengths = [size(bd.fields,1) for bd in [data1, data2]]
            inds = findall(field_lengths .> 0)
            field_lengths = field_lengths[inds]

            # need to use the smallest field length to calculate the fourier_modes of both boundarybasis'
            field_length = minimum(field_lengths)
            basis_order = Int(floor((field_length - 1)/2 ))

        else basis_order = basislength_to_basisorder(PhysicalMedium{2,1},mode_lengths[1])
        end

        @reset sim.method.modal_method.basis_order = basis_order
    end


    # We need the fourier_modes for the boundarybasis
    basis_vec = map(sim.boundarybasis1.basis) do b
        if isempty(b.fourier_modes)
            isempty(b.fields) ? b : fields_to_fouriermodes(b,basis_order)
        else
            b
        end
    end
    @reset sim.boundarybasis1 = BoundaryBasis(basis_vec)

    basis_vec = map(sim.boundarybasis2.basis) do b
        if isempty(b.fourier_modes)
            isempty(b.fields) ? b : fields_to_fouriermodes(b,basis_order)
        else
            b
        end
    end
    @reset sim.boundarybasis2 = BoundaryBasis(basis_vec)

    # We need the fourier_modes for boundarydata if a basis is not provided
    if isempty(sim.boundarybasis1)

        if isempty(sim.boundarydata1.fourier_modes)
            println("The Fourier modes for boundarydata1 are empty, they will be calculated from the fields provided")

            @reset sim.boundarydata1 = fields_to_fouriermodes(sim.boundarydata1, basis_order)
        end    
    end

    if isempty(sim.boundarybasis2)
        if isempty(sim.boundarydata2.fourier_modes)
            println("The Fourier modes for boundarydata2 are empty, they will be calculated from the fields provided")

            @reset sim.boundarydata2 = fields_to_fouriermodes(sim.boundarydata2, basis_order)
        end
    end    

    return sim
end

function point_contact_boundary_data( θs::Vector{T}, bearing::RollerBearing, bc ::BoundaryCondition{TractionType}
    , basis_order ::Int) where{T}

    Z =bearing.number_of_rollers
    μ=bearing.medium.friction_coefficient
    n_order=basis_order

    size=2*n_order*Z+1 

    fourier_coef_p=zeros(size)

    for n in -n_order:n_order
         fourier_coef_p[Int(n_order*Z)+1+Int(n*Z)]=Z/2pi
    end

    fourier_coef_s=μ.*fourier_coef_p 
    bd =  BoundaryData(bc, θs=θs, fourier_modes=hcat(fourier_coef_p,fourier_coef_s))

    return bd
    
end


# boundary_data(sim, loading_profile::BoundaryData{BC})

function BoundaryData(loading_profile::BoundaryData{BC}, bearing::RollerBearing,
      basis_order ::Int, frequency_order:: Int) where{T, BC <: BoundaryCondition{TractionType}}

    
    bc=loading_profile.boundarytype
    θs = loading_profile.θs
    load = loading_profile.fields
    
    # if length(θs)!=length(load)
    #     error("the load and θs must have the same size")
    # end

    Z = bearing.number_of_rollers
    
    #basis_order = basislength_to_basisorder(PhysicalMedium{2,1}, size(loading_profile.fourier_modes,1))


    fouriermodesp=fields_to_fouriermodes(θs,load[:,1], basis_order)
    fouriermodess=fields_to_fouriermodes(θs,load[:,2], basis_order)
    
    
    fouriermodesp=(Z/(2pi)).*fouriermodesp
    fouriermodess=(Z/(2pi)).*fouriermodess
    m=frequency_order

    

    size=2*(basis_order+Z*m)+1 

    fourier_coef_p=zeros(ComplexF64,size) 
    fourier_coef_s=zeros(ComplexF64,size) 

    fourier_coef_p[2*Z*m+1:end]= fouriermodesp
    fourier_coef_s[2*Z*m+1:end]= fouriermodess
    
    bd =  BoundaryData(bc, θs=θs, fourier_modes=hcat(fourier_coef_p,fourier_coef_s))

    return bd
    

end