# module ElasticWaves

# # wave_types.jl
# export field, fields_to_fouriermodes

# # utils.jl and signal_processing.jl
# export estimate_basisorder, fields_to_fouriermodes

using MultipleScattering
using SpecialFunctions
using Statistics
using LinearAlgebra
using StaticArrays: SVector

include("wave_types.jl")
include("bearing.jl")
include("fields.jl")
include("signal_processing.jl")

include("cylindrical/utils.jl")
include("cylindrical/boundary_conditions.jl")
include("cylindrical/fields.jl")

# end # module
