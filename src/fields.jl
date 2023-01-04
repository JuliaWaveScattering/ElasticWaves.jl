import MultipleScattering: field
import MultipleScattering: Shape

"""
    field(potential::HelmholtzPotential, x::AbstractVector)

Returns the value of the potentical at the point `x`
"""
function field(potential::HelmholtzPotential{Dim}, x::AbstractVector{T}) where {Dim, T<:AbstractFloat}

    k = potential.wavenumber
    r, θ = cartesian_to_radial_coordinates(x)

    coefs = permutedims(potential.coefficients, (2,1))

    ms = -potential.basis_order:potential.basis_order
    exps = exp.(im * θ .* ms)

    j = sum(coefs[:,1] .* besselj.(ms,k*r) .* exps)
    h = sum(coefs[:,2] .* hankelh1.(ms,k*r) .* exps)

   return j + h
end

function field(potential::HelmholtzPotential, bearing::RollerBearing; kws...)

    inner_circle = Circle(bearing.inner_radius)
    outer_circle = Circle(bearing.outer_radius)

    return field(potential, outer_circle;
        exclude_region = inner_circle,
        kws...
    )
end

function field(potential::HelmholtzPotential, sh::Shape; kws...)

    x_vec, inds = points_in_shape(sh; kws...)
    # x_vec is a square grid of points and x_vec[inds] are the points in the region.

    xs = x_vec[inds]
    field_mat = zeros(Complex{Float64},length(x_vec))

    fs = [field(potential,x) for x in xs];
    field_mat[inds] = fs

    return  FrequencySimulationResult(field_mat, x_vec, [ω])
end
