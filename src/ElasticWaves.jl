module ElasticWaves

# wave_types.jl
export Elasticity, HelmholtzPotential, DisplacementType, TractionType
export ElasticWave

# bearing.jl
export RollerBearing, BoundaryCondition, DisplacementBoundary, TractionBoundary
export BoundaryData, BoundaryBasis, BearingSimulation

# fields.jl
export field

# cylindrical/fields.jl
export field_modes, displacement, traction, pressure_field_mode, shear_field_mode

# signal_processing.jl
export fouriermodes_to_fields, fields_to_fouriermodes

# cylindrical/utils.jl
export estimate_basisorder # soon to be removed

# cylindrical/elastic_wave.jl
export boundarycondition_mode, boundarycondition_system # soon to be removed
export ModalMethod, PriorMethod, GapMethod


using MultipleScattering
using SpecialFunctions
using LinearAlgebra
using StaticArrays: SVector

include("wave_types.jl")
include("bearing.jl")
include("fields.jl")
include("signal_processing.jl")

include("cylindrical/utils.jl")
include("cylindrical/elastic_wave.jl")
include("cylindrical/fields.jl")

end # module
