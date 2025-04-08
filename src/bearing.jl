struct RollerBearing{T <: AbstractFloat}
    medium::Elastic{2,T} # defining medium
    inner_radius::Union{T,Complex{T}} # due to nondimensionalisation radius might be complex
    "vector of angles delimiting gaps in the inner radius"
    inner_gaps::Vector{T}
    outer_radius::Union{T,Complex{T}}
    "vector of angles delimiting gaps in the outer radius"
    outer_gaps::Vector{T}
    number_of_rollers::Int
    rollers_inside::Bool
    roller_radius::Float64
    "the circumferential distance between the boundaries of adjacent rollers."
    roller_separation::Float64
    "Assumes that the traction on the surface decays like (1/σ) * exp(-π θ^2 / σ^2). The variabel below is σ."
    roller_contact_angular_spread::Float64
    angular_speed::Float64
end

function RollerBearing(; medium::Elastic{2,T},
        inner_radius::Union{T,Complex{T}} = 0.0, 
        outer_radius::Union{T,Complex{T}} = 0.0,
        inner_gaps::Vector{T} = typeof(medium).parameters[2][],
        outer_gaps::Vector{T} = typeof(medium).parameters[2][],
        number_of_rollers::Int = 8,
        roller_radius::T = inner_radius / (typeof(inner_radius)(1.0) + 2.5*number_of_rollers / (2pi)),
        rollers_inside::Bool = true,
        roller_separation::T = 2pi * (rollers_inside ? (inner_radius - roller_radius) : (outer_radius + roller_radius)) / number_of_rollers - 2 * roller_radius,
        roller_contact_angular_spread = 2pi*inner_radius / (5.0 * number_of_rollers), 
        angular_speed::T = 1.0
    ) where T <: AbstractFloat

    if isodd(length(inner_gaps)) && isodd(length(outer_gaps))
        @error "both inner_gaps and outer_gaps need to be an even number of angles"
    end

    if roller_separation < 0
        @error "The roller_separation can not be negative."
    end    

    if rollers_inside && number_of_rollers * (2*roller_radius + roller_separation) > 2pi * (inner_radius - roller_radius) 
        @error "The rollers do not fit! Need to decrease the size of roller_radius or roller_separation or number_of_rollers"
    end

    if !rollers_inside && number_of_rollers * (2*roller_radius + roller_separation) > 2pi * (outer_radius + roller_radius) 
        @error "The rollers do not fit! Need to decrease the size of roller_radius or roller_separation or number_of_rollers"
    end



    println("Note that the traction on the inner boundary (τn,τt) is interpreted to be the vector τ = - τn * er - τt * eθ.")

    # we assume the rollers are in contact with the raceway, but have a neglible indentation into the raceway. In which case we should have that 
    
    # the radius of the circle passing through the centre of the rollers
    R = rollers_inside ? (inner_radius - roller_radius) : (outer_radius + roller_radius)
    d = 2 * roller_radius

    if number_of_rollers >0 && !(2pi * R ≈ number_of_rollers * (d + roller_separation))
        error(" the number of rollers, and separation distance, do not match the geometry of the raceway.")
    end    

    RollerBearing{T}(medium, inner_radius, inner_gaps, outer_radius, outer_gaps, number_of_rollers,rollers_inside, roller_radius, roller_separation, roller_contact_angular_spread, angular_speed)
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

Gives all the information on one boundery of a `RollerBearing`. You can specify the forcing either be giving the `coefficients`, or by specifying the `fields`. The `fields` represent either the displacement or traction when varying the polar angle `θs`.

The relationship between the `fields` and the `coefficients` is given by ``f_{j}(\\theta) = \\sum_m \\text{mode}_{mj} e^{i \\theta m}``, where ``f`` represents the `fields`. To calculate the `fields` from the `coefficients` we can use:
```julia
exps = [
    exp(im * θ * m)
for θ in θs, m = modes];

fields = exps * coefficients
```

If the `coefficients` were specified, then they will be used directly. If the `coefficients` were not given, then `fields` will be used to calculate the `coefficients`: `coefficients = exps \\ fields`.
"""
struct BoundaryData{BC <: BoundaryCondition, T} <: AbstractBoundaryData{BC}
    boundarytype::BC
    θs::Vector{T}
    fields::Matrix{Complex{T}}
    coefficients::Matrix{Complex{T}}
    modes::Vector{Int}

    function BoundaryData(boundarytype::BC, θs::AbstractVector{T}, fields::Matrix{Complex{T}}, coefficients::Matrix{Complex{T}}, modes::AbstractVector{Int}) where {BC <: BoundaryCondition, T}

        if !isempty(coefficients)
            if size(coefficients,2) != 2
                @error "the coefficients should have only two columns" 
            end
    
            is = sortperm_modes(modes);
            if modes != modes[is]
                @warn "A non standard order for the modes has been given. The standard order is assumed to be true by many functions"
            end    
        end    
    
        if !isempty(fields) 
            if size(fields,2) != 2
                @error "the fields should have only two columns"
            end
        end

        return new{BC, T}(boundarytype, θs |> collect, Complex.(fields), coefficients, modes |> collect)
    end

end

function BoundaryData(boundarytype::BC;
        θs::AbstractVector{T} = Float64[],
        fields::Matrix = reshape(Complex{Float64}[],0,2),
        coefficients::Matrix = reshape(Complex{Float64}[],0,2),
        modes::AbstractVector{Int} = Int[]
    ) where {BC <: BoundaryCondition, T}

    if !isempty(coefficients)
        if size(coefficients,2) != 2
            @error "the coefficients should have only two columns" 
        end

        if isempty(modes)
            order = basislength_to_basisorder(PhysicalMedium{2,1},size(coefficients,1))
            modes = -order:order

        elseif size(coefficients,1) != length(modes)
            @error "the Fourier coefficients should have the same length as modes"
        end

        is = sortperm_modes(modes);
        modes = modes[is];
        coefficients = coefficients[is,:]; 
    end   

    if !isempty(fields)
        
        if size(fields,2) != 2
            @error "the fields should have only two columns"
        end
        if size(fields,1) != length(θs)
            @error "the fields should have the same length as θs"
        end

        is = sortperm(θs);
        θs = θs[is];
        fields = fields[is,:];
    end

    return BoundaryData(boundarytype, θs |> collect, Complex.(fields), Complex.(coefficients), modes |> collect)
end

"""
    select_modes(bd::Union{BoundaryData,HelmholtzPotential}, modes::AbstractVector{Int})

Returns a 'BoundaryData' or 'HelmholtzPotential' with only the 'modes' given.    
"""
function select_modes(bd::Union{BD,H}, modes::AbstractVector{Int}) where {BD <: BoundaryData, H <: HelmholtzPotential}
    
    if !(setdiff(modes,bd.modes) |> isempty)
        error("There are not enough modes in the boundarydata given to select the modes = $modes")
    end
    
    inds = [findfirst(m .== bd.modes) for m in modes]
    
    @reset bd.modes = modes
    
    if typeof(bd) <: HelmholtzPotential
        @reset bd.coefficients = bd.coefficients[:,inds]
    else    
        @reset bd.coefficients = bd.coefficients[inds,:]
    end    

    return bd
end

"""
    select_modes(wave::ElasticWave, modes::AbstractVector{Int})

Returns an 'ElasticWave' with only the 'modes' given.    
"""
function select_modes(wave::ElasticWave, modes::AbstractVector{Int})
    
    @reset wave.potentials = [select_modes(p, modes) for p in wave.potentials]
    return wave
end

"""
    BoundaryData(wave::ElasticWave{2}, radius, wave::ElasticWave{2})

A convenience function to predict the boundary data  according to the solution 'wave' using the modes only.    
"""
function BoundaryData(boundarycondition::BoundaryCondition, radius::AbstractFloat, wave::ElasticWave{2})

    fs = field_modes(wave, radius, boundarycondition.fieldtype)

    return BoundaryData(boundarycondition; coefficients = fs, modes = wave.potentials[1].modes)
end

"""
    BoundaryData(boundarydata::BoundaryData, radius, wave::ElasticWave{2})

A convenience function to predict the data on a boundary according to the solution 'wave'. The format of the data, including the mesh and modes, is determined from boundarydata.
"""
function BoundaryData(bd::BoundaryData{BC,T}, radius::T, wave::ElasticWave{2}) where {T <: AbstractFloat, BC}
    
    θs = bd.θs
    modes = bd.modes
    
    # Calculate the fields according to wave
    xs = [
        radial_to_cartesian_coordinates([radius,θ]) 
    for θ in θs];
    
    fs = [field(wave, x, bd.boundarytype.fieldtype) for x in xs];
    fs = hcat(fs...) |> transpose |> collect
    fs = reshape(fs,:,2)

    # Calculate the fourier modes according to wave
    modes_setdiff = setdiff(modes, wave.potentials[1].modes);
    if modes_setdiff |> length > 0
        @warn "The modes in wave are not sufficient to calculate all the modes from the BoundaryData. Will set coefficients of the missing to zero."
    end

    # using only the modes which are shared, the others are zero
    modes_intersect = intersect(modes, wave.potentials[1].modes);
    wave_intersect = select_modes(wave,modes_intersect)

    f_modes = field_modes(wave_intersect, radius, bd.boundarytype.fieldtype)
    f_modes = reshape(f_modes,:,2)

    modes_all = [modes_setdiff; modes_intersect]
    f_modes = vcat(zeros(Complex{T},length(modes_setdiff),2), f_modes)
   
    inds = sortperm_modes(modes_all)
    
    return BoundaryData(bd.boundarytype; 
        θs = θs, fields = Matrix{Complex{T}}(fs), 
        modes = modes_all[inds], coefficients = Matrix{Complex{T}}(f_modes[inds,:])
    )
end

struct BoundaryBasis{BD <: AbstractBoundaryData}
    basis::Vector{BD}

    function BoundaryBasis(basis::Vector{BD}) where {BD <: AbstractBoundaryData}

        T = basis[1].coefficients |> typeof
        Tc = T.parameters[1]
        
        T = basis[1].fields |> typeof
        Tf = T.parameters[1]
    
        # find all the angles and modes represented. 
        θs = union([b.θs for b in basis]...)
        modes = union([b.modes for b in basis]...)
    
        # Pad each element of the basis to cover the same angles and modes
        basis_padded = map(basis) do b
            add_modes = setdiff(modes, b.modes)
            b_modes = [add_modes; b.modes]
            b_coes = [zeros(Tc, add_modes |> length , 2); b.coefficients]
    
            inds = sortperm_modes(b_modes)
            @reset b.modes = b_modes[inds]
            @reset b.coefficients = b_coes[inds,:]
    
            add_θs = setdiff(θs, b.θs)
            b_θs = [add_θs; b.θs]
            b_fields = [zeros(Tf, add_θs |> length , 2); b.fields]
    
            inds = sortperm(b_θs)
            @reset b.θs = b_θs[inds]
            @reset b.fields = b_fields[inds,:]
    
            b
        end    
    
        return new{BD}(basis_padded)
    end
end

import Base: isempty
isempty(bd::BoundaryData) = all([isempty(bd.fields),isempty(bd.coefficients)])  
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

mutable struct BearingSimulation{M <: SolutionMethod, BC1 <: BoundaryCondition, BC2 <: BoundaryCondition, BD1 <: BoundaryData, BD2 <: BoundaryData, T}
    ω::T
    method::M
    nondimensionalise::Bool 
    bearing::RollerBearing{T}
    boundarydata1::BoundaryData{BC1,T}
    boundarydata2::BoundaryData{BC2,T}
    boundarybasis1::BoundaryBasis{BD1}
    boundarybasis2::BoundaryBasis{BD2}
end

function BearingSimulation(ω::Real, method::M, bearing::RollerBearing{T}, 
        boundarydata1::BoundaryData{BC1,T},
        boundarydata2::BoundaryData{BC2,T};
        nondimensionalise::Bool = true,
        boundarybasis1::BoundaryBasis{BD1}  = BoundaryBasis([BoundaryData(boundarydata1.boundarytype)]),
        boundarybasis2::BoundaryBasis{BD2}  = BoundaryBasis([BoundaryData(boundarydata2.boundarytype)])
    ) where {T, M <: SolutionMethod, BC1 <: BoundaryCondition, BC2 <: BoundaryCondition, BD1 <: BoundaryData, BD2 <: BoundaryData,}

    sim = BearingSimulation{M,BC1,BC2,BD1,BD2,T}(T(ω), method, nondimensionalise, bearing, boundarydata1, boundarydata2, boundarybasis1, boundarybasis2)

    return setup(sim)
end

import MultipleScattering: boundary_data
"""
    boundary_data(wave::ElasticWave{2}, sim::BearingSimulation)

A convenience function to calculate the boundary data of 'sim' according to the solution 'wave'.    
"""
function boundary_data(sim::BearingSimulation, wave::ElasticWave{2})

    if wave.potentials[1].modes |> isempty 
        @error "The wave potentials are empty. Can not calculate field."
    end     
    bd1 = sim.boundarydata1;
    bd2 = sim.boundarydata2;

    radius_inner_outer(in) = (in == true) ? sim.bearing.inner_radius : sim.bearing.outer_radius

    inners = [bd1.boundarytype.inner, bd2.boundarytype.inner]
    radiuses = radius_inner_outer.(inners)
    
    
    bd1_location = BoundaryData(bd1.boundarytype; θs = bd1.θs, modes = bd1.modes)
    bd2_location = BoundaryData(bd2.boundarytype; θs = bd2.θs, modes = bd2.modes)

    bd1 = BoundaryData(bd1_location, radiuses[1], wave)
    bd2 = BoundaryData(bd2_location, radiuses[2], wave)

    return [bd1, bd2]
end

function nondimensionalise!(boundarydata::BoundaryData{BoundaryCondition{TractionType}}, sim::BearingSimulation)

    λ2μ = sim.bearing.medium.ρ * sim.bearing.medium.cp^2
    boundarydata.fields[:] = boundarydata.fields ./ λ2μ
    boundarydata.coefficients[:] = boundarydata.coefficients ./ λ2μ

    return boundarydata
end

function nondimensionalise!(boundarydata::BoundaryData{BoundaryCondition{DisplacementType}}, sim::BearingSimulation)

    kp = sim.ω / sim.bearing.medium.cp;
    boundarydata.fields[:] = boundarydata.fields .* kp
    boundarydata.coefficients[:] = boundarydata.coefficients .* kp

    return boundarydata
end

function nondimensionalise(bearing::RollerBearing, ω::Number)

    kp = ω / bearing.medium.cp;
    non_medium = Elastic(2; ρ = 1 / ω^2, cp = ω, cs = bearing.medium.cs * kp)
    
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
    
    # just normalise each of elements of the basis as their size does not matter
    # normalize!(sim.boundarybasis1)
    # normalize!(sim.boundarybasis2)

    # nondimensionalise bearing
    sim.bearing = nondimensionalise(sim.bearing,sim.ω)

    return sim
end


function BearingSimulation(ω::Real, bearing::RollerBearing{T}, boundarydata1::BD1, boundarydata2::BD2;
        method::M = NoBearingMethod(), kws...
    ) where {T, M <: SolutionMethod, BD1 <: BoundaryData, BD2 <: BoundaryData}

    if typeof(method) == NoBearingMethod
        @warn "no method was chosen for this bearing simulation. Will choose the default ModalMethod which makes no assumption about the forcing or bearings"
        
        method = ModalMethod()
    end    
    
    return BearingSimulation(ω, method, bearing, boundarydata1, boundarydata2; kws...) 
end

function setup(sim::BearingSimulation{ModalMethod})
    
    max_mode = sim.method.maximum_mode
    modes = sim.method.modes
    boundarydata1 = sim.boundarydata1
    boundarydata2 = sim.boundarydata2

    get_order(coefficients) = if isempty(coefficients)
        - 1 
    else 
        l = size(coefficients,1)
        l = 2 * Int(floor((l - 1) / 2)) + 1
        basislength_to_basisorder(PhysicalMedium{2,1}, l)
    end    

    # if the modes not given, then determine the modes from the data
    if isempty(modes)

        bd_to_modes(bd) = isempty(bd.modes) ? (-get_order(bd.fields):get_order(bd.fields)) : bd.modes

        modes_vec = bd_to_modes.([boundarydata1, boundarydata2])

        modes = intersect(modes_vec...) |> collect

        if isempty(modes)
            @error "boundarydata1 and boundarydata2 share no modes in common. Fourier modes are required for the ModalMethod"
        end

        is = findall(abs.(modes) .< max_mode)
        modes = modes[is]

        is = sortperm_modes(modes);
        modes = modes[is];

        @reset sim.method.modes = modes
    end

   # we need the fourier modes of the boundary data
    if isempty(boundarydata1.coefficients) 
        # || (get_order(boundarydata1.fields) > get_order(boundarydata1.coefficients) && get_order(boundarydata1.coefficients) < basis_order )
        # println("The Fourier modes for boundarydata1 are empty, or the fields have more data. The Fourier modes will be calculated from the fields provided.")
        println("The Fourier coefficients for boundarydata1 are empty. The Fourier coefficients will be calculated from the fields provided.")
        boundarydata1 = fields_to_fouriermodes(boundarydata1, modes)
        
    elseif modes != boundarydata1.modes 
        boundarydata1 = select_modes(boundarydata1, modes)
    end

    if isempty(boundarydata2.coefficients) 
        println("The Fourier coefficients for boundarydata2 are empty. The Fourier coefficients will be calculated from the fields provided.")
        boundarydata2 = fields_to_fouriermodes(boundarydata2, modes)
    
    elseif modes != boundarydata2.modes 
        boundarydata2 = select_modes(boundarydata2, modes)
    end

    @reset sim.boundarydata1 = boundarydata1
    @reset sim.boundarydata2 = boundarydata2

    return sim
end

function setup(sim::BearingSimulation{PriorMethod})

    modes = sim.method.modal_method.modes;

    if isempty(sim.boundarybasis1) && isempty(sim.boundarybasis2)
        error("The PriorMethod requires a boundary basis either for the first boundary (keyword 'boundary_basis1') or the second boundary (keyword 'boundary_basis2') ")
    end    

    # we need the fields of the boundary data for the PriorMethod
    if isempty(sim.boundarydata1.fields) && isempty(sim.boundarydata2.fields)

        error("The fields of boundarydata1 or boundarydata2 are empty. We expect that the fields are given for the PriorMethod because as the number of coefficients can be far greater than the length of the fields. Meaning it is not possible to calculate the coefficients of the boundary data directly")
    end

    if modes |> isempty
        println("No modes were specified. Will use the as many modes as possible according to the data provided (boundarybasis1, boundarybasis2, boundarydata1, and boundarydata2)")

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

        basis_boundary_number = isempty(sim.boundarybasis1) ? 2 : 1

        get_order(coefficients) = isempty(coefficients) ? - 1 : basislength_to_basisorder(PhysicalMedium{2,1}, size(coefficients,1))
        bd_to_modes(bd) = isempty(bd.modes) ? (-get_order(bd.fields):get_order(bd.fields)) : bd.modes

        modes_vec = bd_to_modes.([data1, data2])
        modes = intersect(modes_vec...) |> collect

        if length(intersect(modes,modes_vec[basis_boundary_number])) < length(modes_vec[basis_boundary_number])
            error("There is not enough boundary data for the modes specified in the boundarybasis. Try increasing the amount of data on the boundarydata.") 
        end

        # use the standard ordering
        modes = modes[sortperm_modes(modes)]

        @reset sim.method.modal_method.modes = modes
    end

    # We need the coefficients for the boundarybasis
    basis_vec = map(sim.boundarybasis1.basis) do b
        if isempty(b.coefficients)
            isempty(b.fields) ? b : fields_to_fouriermodes(b,modes)
        else
            select_modes(b, modes)
        end
    end
    @reset sim.boundarybasis1 = BoundaryBasis(basis_vec)

    basis_vec = map(sim.boundarybasis2.basis) do b
        if isempty(b.coefficients)
            isempty(b.fields) ? b : fields_to_fouriermodes(b,modes)
        else
            select_modes(b, modes)
        end
    end
    @reset sim.boundarybasis2 = BoundaryBasis(basis_vec)

    # We need the coefficients for boundarydata if a basis is not provided
    if isempty(sim.boundarybasis1)

        if isempty(sim.boundarydata1.coefficients)
            println("The Fourier modes for boundarydata1 are empty, they will be calculated from the fields provided")

            @reset sim.boundarydata1 = fields_to_fouriermodes(sim.boundarydata1, modes)
        end
    end

    if isempty(sim.boundarybasis2)
        if isempty(sim.boundarydata2.coefficients)
            println("The Fourier modes for boundarydata2 are empty, they will be calculated from the fields provided")

            @reset sim.boundarydata2 = fields_to_fouriermodes(sim.boundarydata2, modes)
        end
    end    

    return sim
end
