module ElasticWaves

# elasticity.jl
export Elastic, ElasticWave
export regular_basis_function, HelmholtzPotential 
export DisplacementType, TractionType

# source.jl
export plane_z_shear_source

# bearing.jl
export RollerBearing, BoundaryCondition, DisplacementBoundary, TractionBoundary
export BoundaryData, BoundaryBasis, BearingSimulation

# fields.jl
export field

# spherical/fields.jl
export pressure_regular_basis, shearΦ_regular_basis, shearχ_regular_basis
export pressure_outgoing_basis, shearΦ_outgoing_basis, shearχ_outgoing_basis

# spherical/t-matrix.jl
export modal_system, t_matrix, internal_matrix

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
using BlockDiagonals


include("elasticity.jl")
include("source.jl")

include("bearing.jl")
include("fields.jl")
include("signal_processing.jl")

include("spherical/fields.jl")
include("spherical/t-matrix.jl")

include("cylindrical/utils.jl")
include("cylindrical/elastic_wave.jl")
include("cylindrical/fields.jl")

end # module
