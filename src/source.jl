# Create elastic wave sources

"""
    plane_z_source(medium::Elastic{3,T}, pos::AbstractArray{T}, amplitude::Union{T,Complex{T}}

Creates an incident plane wave for the displacement in the form ``A e^{i k z} (0,1,0)``. The coefficients of the Debye potentials which correspond to this plane wave are given in "Resonance theory of elastic waves ultrasonically scattered from an elastic sphere - 1987".
"""
function plane_z_shear_source(medium::Elastic{3,T}, pos::AbstractArray{T} = zeros(T,3),
            amplitude::Union{T,Complex{T}} = one(T)
        ) where {T}

    S = PlanarSymmetry{3}

    # code assumes wave propagates in z-direction and polarised in y-direction     
    direction = [zero(T),zero(T),one(T)]
    direction = direction / norm(direction)

    polarisation = [zero(T),one(T),zero(T)]
    polarisation = polarisation / norm(polarisation)

    # Convert to SVector for efficiency and consistency
    # position = SVector(position...)

    function source_field(x,ω)
        # x_width = norm((x - position) - dot(x - position, direction)*direction)
        
        return amplitude .* exp(im * ω / medium.cs * dot(x - pos, direction)) .* polarisation
    end

    function spherical_expansion(order,centre,ω)
        
        ks2 = (ω / medium.cs)^2

        pcoefs = [Complex{T}(0.0) for l = 0:order for m = -l:l] |> transpose
        p_potential = HelmholtzPotential{3}(medium.cp, ω / medium.cp, [pcoefs; 0 .* pcoefs])
        
        Φcoefs = - T(sqrt(pi)) * sum(source_field(centre,ω) .* polarisation) .*
        [
            (m == 1 | m == -1) ?  Complex{T}(m * (-1.0im)^l * sqrt(T(2l + 1) / T((1+l)*l))) : Complex{T}(0)
        for l = 0:order for m = -l:l] |> transpose
        Φ_potential = HelmholtzPotential{3}(medium.cs, ω / medium.cs, [Φcoefs; 0 .* Φcoefs])
        
        χcoefs = - Complex{T}(sqrt(pi) / ks2) * sum(source_field(centre,ω) .* polarisation) .*
        [
            (m == 1 | m == -1) ?  Complex{T}((-1.0im)^l * sqrt(T(2l + 1) / T((1+l)*l))) : Complex{T}(0)
        for l = 0:order for m = -l:l] |> transpose
        χ_potential = HelmholtzPotential{3}(medium.cs, ω / medium.cs, [χcoefs; 0 .* χcoefs])


        return ElasticWave(ω, medium, [p_potential,Φ_potential,χ_potential])
    end

    return RegularSource{Elastic{3,T},S}(medium, source_field, spherical_expansion)
end

# vs = regular_basis_function(source.medium, ω)
# regular_coefficients = regular_spherical_coefficients(source)

# for x close to centre
# source_field(x,ω) ~= sum(regular_coefficients(basis_order,centre,ω) .* vs(basis_order, x - centre))
