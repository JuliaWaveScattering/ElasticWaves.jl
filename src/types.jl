
## field types
"""
    FieldType

A type used to specify what type of physical field, such as traction or displacement.
"""
abstract type FieldType end

struct DisplacementType <: FieldType end
struct TractionType <: FieldType end

## methods to solve for waves in bearings 

abstract type SolutionMethod end

abstract type BearingMethod <: SolutionMethod end

struct NoBearingMethod <: SolutionMethod end

struct ModalMethod <: SolutionMethod
    tol::Float64
    # to use Tikhonov regularization give a non-zero parameter
    regularisation_parameter::Float64
    only_stable_modes::Bool
    maximum_mode::Int
    modes::Vector{Int}
    mode_errors::Vector{Float64}

    function ModalMethod(tol::Float64, regularisation_parameter::Float64, only_stable_modes::Bool, maximum_mode::Int, modes::Vector{Int}, mode_errors::Vector{Float64})

        is = findall(abs.(modes) .< maximum_mode)
        modes = modes[is]

        is = sortperm_modes(modes);
        modes = modes[is]

        if !isempty(mode_errors)
            if length(mode_errors) != length(modes)
                @warn("both mode_errors and modes where given but have different lengths.")

            else mode_errors = mode_errors[is]
            end
        end

        new(tol, regularisation_parameter, only_stable_modes, maximum_mode, modes, mode_errors)
    end    
end

struct GapMethod <: BearingMethod end

struct PriorMethod <: SolutionMethod
    tol::Float64
    # to use Tikhonov regularization give a non-zero parameter
    regularisation_parameter::Float64
    modal_method::ModalMethod
    condition_number::Float64
    boundary_error::Float64
end

function ModalMethod(; 
        tol::Float64 = eps(Float64)^(1/2), 
        regularisation_parameter::Float64 = eps(typeof(tol)) / tol,
        only_stable_modes::Bool = true,
        maximum_mode::Int = 160,
        modes::Vector{Int} = Int[],
        mode_errors::Vector = Float64[]
    )

    if !only_stable_modes 
        @warn "only_stable_modes was set to false. This means that potentially ill-posed (or unstable) modes will attempt to be solved, which could lead to non-sense solutions." 
    end

    ModalMethod(tol, regularisation_parameter, only_stable_modes, maximum_mode, modes, mode_errors)
end

function PriorMethod(; 
        tol::Float64 = eps(Float64)^(1/2), 
        regularisation_parameter::Float64 = zero(Float64),
        modes::Vector{Int} = Int[],
        modal_method = ModalMethod(tol = tol, modes = modes),
        condition_number = -one(Float64),
        boundary_error = -one(Float64),

    )
    PriorMethod(tol, regularisation_parameter, modal_method, condition_number, boundary_error)
end