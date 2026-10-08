
## field types
"""
    FieldType

A type used to specify what type of physical field, such as traction or displacement.
"""
abstract type FieldType end

"""
    DisplacementType <: FieldType

Used to ask for, or to specify, the displacement field. See also [`TractionType`](@ref).
"""
struct DisplacementType <: FieldType end

"""
    TractionType <: FieldType

Used to ask for, or to specify, the traction field. See also [`DisplacementType`](@ref).
"""
struct TractionType <: FieldType end

"""
    PotentialType <: FieldType

Used to work directly with the Helmholtz potentials, without specifying a physical field.
"""
struct PotentialType <: FieldType end

## methods to solve for waves in bearings 

abstract type SolutionMethod end
abstract type BearingMethod <: SolutionMethod end
abstract type AbstractPriorMethod <: SolutionMethod end

struct NoBearingMethod <: SolutionMethod end

"""
    ModalMethod(; tol = sqrt(eps(Float64)), only_stable_modes = true, maximum_mode = 160, modes = Int[])

Solve for the waves in a [`RollerBearing`](@ref) one Fourier mode at a time. For each mode the two boundary conditions of a [`BearingSimulation`](@ref) lead to a 4 x 4 linear system for the coefficients of the pressure and shear potentials.

A mode is considered unstable when its estimated relative error is larger than `tol`. If `only_stable_modes = true`, the modes are solved in order of increasing `abs(mode)` until the first unstable mode is found, and then that mode and all higher modes are discarded. If `only_stable_modes = false`, all modes are kept, which can lead to very large errors.

If `modes` is empty, the modes are determined from the boundary data, keeping only `abs(mode) < maximum_mode`.

After solving with [`ElasticWave(::BearingSimulation)`](@ref), the method stored in the wave holds the `modes` that were kept and the estimated error of each in `mode_errors`.
"""
struct ModalMethod <: SolutionMethod
    tol::Float64
    only_stable_modes::Bool
    maximum_mode::Int
    modes::Vector{Int}
    mode_errors::Vector{Float64}

    function ModalMethod(tol::Float64, only_stable_modes::Bool, maximum_mode::Int, modes::AbstractVector{Int}, mode_errors::Vector{Float64})

        is = findall(abs.(modes) .< maximum_mode)
        modes = modes[is] |> collect

        # is = sortperm_modes(modes);
        # modes = modes[is]

        if !isempty(mode_errors)
            if length(mode_errors) != length(modes)
                @warn("both mode_errors and modes where given but have different lengths.")

            else mode_errors = mode_errors[is]
            end
        end

        new(tol, only_stable_modes, maximum_mode, modes, mode_errors)
    end    
end

struct GapMethod <: BearingMethod end

"""
    PriorMethod(; tol = sqrt(eps(Float64)), modes = Int[], modal_method = ModalMethod(tol = tol, modes = modes))

Solve for the waves in a [`RollerBearing`](@ref) by using prior knowledge about the boundary data. The prior is given as a [`BoundaryBasis`](@ref) when defining the [`BearingSimulation`](@ref): the data on that boundary is assumed to be a linear combination of the elements of the basis, and the coefficients of this combination are fitted to the `fields` of the measured boundary data.

This is typically used to solve inverse problems where there are only a few measurements. For example, to predict the traction on the inner boundary from a few sensors placed on the outer boundary.

The Fourier modes used, and the tolerance for which of them are stable, are given by `modal_method`, see [`ModalMethod`](@ref).

After solving with [`ElasticWave(::BearingSimulation)`](@ref), the method stored in the wave holds the `condition_number` of the system which was solved, and the relative error `boundary_error` of the fit to the measured data.
"""
struct PriorMethod <: AbstractPriorMethod
    tol::Float64
    modal_method::ModalMethod
    condition_number::Float64
    boundary_error::Float64
end

function ModalMethod(; 
        tol::Float64 = eps(Float64)^(1/2), 
        only_stable_modes::Bool = true,
        maximum_mode::Int = 160,
        modes::AbstractVector{Int} = Int[],
        mode_errors::Vector = Float64[]
    )

    if !only_stable_modes 
        @warn "only_stable_modes was set to false. This means that potentially ill-posed (or unstable) modes will attempt to be solved, which could lead to non-sense solutions." 
    end

    ModalMethod(tol, only_stable_modes, maximum_mode, modes, mode_errors)
end

function PriorMethod(; 
        tol::Float64 = eps(Float64)^(1/2), 
        modes::Vector{Int} = Int[],
        modal_method = ModalMethod(tol = tol, modes = modes),
        condition_number = -one(Float64),
        boundary_error = -one(Float64),

    )
    PriorMethod(tol, modal_method, condition_number, boundary_error)
end

"""
    ConstantRollerSpeedMethod(; loading_modes = Int[], ratio_shear_to_normal = 0.0, modal_method, tol)

A prior method which assumes that the rollers move at a constant angular speed, and that the boundary data is due to an unknown loading profile on the rollers. The Fourier modes of the loading profile to solve for are given by `loading_modes`, and several frequencies are solved together with [`ElasticWaveVector`](@ref). See also [`PriorMethod`](@ref).
"""
struct ConstantRollerSpeedMethod <: AbstractPriorMethod
    tol::Float64
    modal_method::ModalMethod
    condition_number::Float64
    boundary_error::Float64
    loading_modes::Vector{Int}
    loading_coefficients::Vector{Complex{Float64}}
    ratio_shear_to_normal::Float64
end

function ConstantRollerSpeedMethod(;
        tol::Float64 = eps(Float64)^(1/2), 
        modal_method = ModalMethod(tol = tol, modes = modes),
        condition_number = -one(Float64),
        boundary_error = -one(Float64),
        loading_modes::AbstractVector{Int} = Int[],
        loading_coefficients::Vector{Complex{Float64}} = Complex{Float64}[],
        ratio_shear_to_normal::Float64 = zero(Float64)
    )
    ConstantRollerSpeedMethod(tol, modal_method, condition_number, boundary_error, loading_modes |> collect, loading_coefficients, ratio_shear_to_normal)
end