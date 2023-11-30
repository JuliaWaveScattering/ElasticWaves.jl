# Create elastic wave sources

"""
    plane_z_source(medium::Elastic{3,T}, position::AbstractArray{T}, amplitude::Union{T,Complex{T}}

Creates an incident plane wave for the displacement in the form ``A e^{i k z} (0,1,0)``. The coefficients of the Debye potentials which correspond to this plane wave are given in "Resonance theory of elastic waves ultrasonically scattered from an elastic sphere - 1987".
"""
function plane_z_source(medium::Elastic{3,T}, position::AbstractArray{T},
            amplitude::Union{T,Complex{T}} = one(T)
        ) where {T}

    # code assumes wave propagates in z-direction and polarised in y-direction     
    direction = [zero(T),zero(T),one(T)]
    polarisation = [zero(T),one(T),zero(T)]


    # Convert to SVector for efficiency and consistency
    # position = SVector(position...)

    function source_field(x,ω)
        x_width = norm((x - position) - dot(x - position, direction)*direction)
        if (causal && dot(x - position,direction) < 0) || x_width > beam_width/2
            zero(Complex{T})
        else
            amplitude .* exp(im*ω/medium.c*dot(x - position, direction)) .* polarisation
        end
    end

    function source_coef(order,centre,ω)
        # plane-wave expansion for complex vectors
        r, θ, φ  = cartesian_to_radial_coordinates(direction)
        Ys = spherical_harmonics(order, θ, φ)
        lm_to_n = lm_to_spherical_harmonic_index

        return T(4pi) * source_field(centre,ω) .*
        [
            Complex{T}(im)^l * (-one(T))^m * Ys[lm_to_n(l,-m)]
        for l = 0:order for m = -l:l]
    end

    return RegularSource{Acoustic{T,3},S}(medium, source_field, source_coef)
end

# vs = regular_basis_function(source.medium, ω)
# regular_coefficients = regular_spherical_coefficients(source)

# for x close to centre
# source_field(x,ω) ~= sum(regular_coefficients(basis_order,centre,ω) .* vs(basis_order, x - centre))
