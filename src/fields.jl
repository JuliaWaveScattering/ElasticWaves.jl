import MultipleScattering: field
import MultipleScattering: Shape

function displacement(wave::ElasticWave, x::Vector{T}) where T <: AbstractFloat
    field(wave, x, DisplacementType())
end

function traction(wave::ElasticWave, x::Vector{T}) where T <: AbstractFloat
    field(wave, x, TractionType())
end

function field(potential::HelmholtzPotential, sh::Shape; kws...)

    x_vec, inds = points_in_shape(sh; kws...)
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.


    xs = x_vec[inds]
    field_mat = zeros(Complex{Float64},length(x_vec))

    fs = [field(potential,x) for x in xs];
    field_mat[inds] = fs

    ω = potential.wavenumber * potential.wavespeed

    if !(ω ≈ real(ω))
        @warn "the angular frequency potential.wavenumber * potential.wavespeed is expected to be real. Will keep only the real part"
    end

    return FrequencySimulationResult(field_mat, x_vec, [real(ω)])
end
