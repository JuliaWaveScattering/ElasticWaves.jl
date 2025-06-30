# Create elastic wave sources
import MultipleScattering: point_source, AbstractSource

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

    polarisation = [one(T), zero(T),zero(T)]
    polarisation = polarisation / norm(polarisation)

    # Convert to SVector for efficiency and consistency
    # position = SVector(position...)

    function source_field(x,ω)
        # x_width = norm((x - position) - dot(x - position, direction)*direction)
        
        return amplitude .* exp(im * ω / medium.cs * dot(x - pos, direction)) .* polarisation
    end

    function spherical_expansion(order,centre,ω)
        
        ks = ω / medium.cs
        ks2 = ks^2

        pcoefs = [Complex{T}(0.0) for l = 0:order for m = -l:l] 
        # p_potential = HelmholtzPotential(medium.cp, ω / medium.cp, [pcoefs; 0 .* pcoefs])
        
        Φcoefs = T(sqrt(pi)) * sum(source_field(centre,ω) .* polarisation) .*
        [
            (abs(m) == 1) ?  Complex{T}(m * (1.0im)^l * sqrt(T(2l + 1) / T((1+l)*l))) : Complex{T}(0)
        for l = 0:order for m = -l:l] 
        # Φ_potential = HelmholtzPotential(medium.cs, ω / medium.cs, [Φcoefs; 0 .* Φcoefs])
        
        χcoefs = T(sqrt(pi)) * sum(source_field(centre,ω) .* polarisation) .*
        [
            (abs(m) == 1) ?  Complex{T}((1.0im)^l * sqrt(T(2l + 1) / T((1+l)*l))) : Complex{T}(0)
        for l = 0:order for m = -l:l] 
        # χ_potential = HelmholtzPotential(medium.cs, ω / medium.cs, [χcoefs; 0 .* χcoefs])

        # return ElasticWave(ω, medium, [pcoefs,Φcoefs,χcoefs])
        return [pcoefs Φcoefs χcoefs] ./ (-im * ks) |> transpose
    end

    return RegularSource{Elastic{3,T},S}(medium, source_field, spherical_expansion)
end

"""
    pressure_point_source

    NOTE: currently untested
"""
function pressure_point_source(medium::Elastic{2,T}, source_position::AbstractVector, amplitude::Union{T,Complex{T},Function} = one(T))::RegularSource{Elastic{2,T}} where T <: AbstractFloat

    # Convert to SVector for efficiency and consistency
    source_position = SVector{2,T}(source_position)

    if typeof(amplitude) <: Number
        amp(ω) = amplitude
    else
        amp = amplitude
    end
    source_field(x,ω) = (amp(ω)*im)/4 * hankelh1(0,ω/medium.cp * norm(x-source_position))

    function source_coef(order,centre,ω)
        k = ω/medium.cp
        r, θ = cartesian_to_radial_coordinates(centre - source_position)

        # using Graf's addition theorem
        return (amp(ω)*im)/4 * [hankelh1(-n,k*r) * exp(-im*n*θ) for n = -order:order]
    end

    return RegularSource{Elastic{2,T},WithoutSymmetry{2}}(medium, source_field, source_coef)
end

"""
    pressure_point_source

    NOTE: currently untested
"""
function shear_point_source(medium::Elastic{2,T}, source_position::AbstractVector, amplitude::Union{T,Complex{T},Function} = one(T))::RegularSource{Elastic{2,T}} where T <: AbstractFloat

    # Convert to SVector for efficiency and consistency
    source_position = SVector{2,T}(source_position)

    if typeof(amplitude) <: Number
        amp(ω) = amplitude
    else
        amp = amplitude
    end
    source_field(x,ω) = (amp(ω)*im)/4 * hankelh1(0,ω/medium.cs * norm(x-source_position))

    function source_coef(order,centre,ω)
        k = ω/medium.cs
        r, θ = cartesian_to_radial_coordinates(centre - source_position)

        # using Graf's addition theorem
        return (amp(ω)*im)/4 * [hankelh1(-n,k*r) * exp(-im*n*θ) for n = -order:order]
    end

    return RegularSource{Elastic{2,T},WithoutSymmetry{2}}(medium, source_field, source_coef)
end
# vs = regular_basis_function(source.medium, ω)
# regular_coefficients = regular_spherical_coefficients(source)

# for x close to centre
# source_field(x,ω) ~= sum(regular_coefficients(basis_order,centre,ω) .* vs(basis_order, x - centre))

struct SourceMap
    locations::AbstractArray{Float64}
    amplitudes::Vector{ComplexF64}
end

# Needs to be checked against a method which directly evalutes how the source contributes to  boundarytype at specific points on the boundary.
function boundary_data(ω::T, bearing::RollerBearing, boundarytype::BC, p_map::SourceMap, s_map::SourceMap, modes::Vector{Int}, θs::AbstractVector = T[]) where {BC <: BoundaryCondition, T}

    order = maximum(modes)
    basis_length = 2*order+1

    locations_p = p_map.locations
    amplitudes_p = p_map.amplitudes
    locations_s = s_map.locations
    amplitudes_s = s_map.amplitudes

    number_of_p_sources = length(amplitudes_p)
    number_of_s_sources = length(amplitudes_s)

    acoustic_medium = Acoustic(spatial_dimension(bearing.medium); ρ = bearing.medium.ρ, c = bearing.medium.cp)
    shear_medium = Acoustic(spatial_dimension(bearing.medium); ρ = bearing.medium.ρ, c = bearing.medium.cs)
    
    # potentials of the form ΣH(kr')J(kr)exp(i*n*(θ-θ'))
    # I checked the expersions below against Jess's
    if boundarytype.inner == true
        p_coes = [
            (-im/4) * amplitudes_p[i] * vcat(
                outgoing_translation_matrix(acoustic_medium, order, 0, ω, -locations_p[i,:]) ,
                zeros(ComplexF64, 1, basis_length)
            )
        for i in 1:number_of_p_sources]
        
        s_coes = [
            (-im/4) * amplitudes_s[i] * vcat(
                outgoing_translation_matrix(shear_medium, order, 0, ω, -locations_s[i,:])  ,
                zeros(ComplexF64, 1, basis_length)
            ) 
        for i in 1:number_of_s_sources]
    else    
        p_coes = [
            (-im/4) * amplitudes_p[i] * vcat(
                zeros(ComplexF64, 1, basis_length),
                regular_translation_matrix(acoustic_medium, order, 0, ω, -locations_p[i,:]) 
            ) 
        for i in 1:number_of_p_sources]
        
        s_coes = [
            (-im/4) * amplitudes_s[i] * vcat(
                zeros(ComplexF64, 1, basis_length),
                regular_translation_matrix(shear_medium, order, 0, ω, -locations_s[i,:])  
            ) 
        for i in 1:number_of_s_sources]
    end

    p_coes = sum(p_coes) 
    s_coes = sum(s_coes) 

    #is=sortperm_modes(modes)

    #modes=modes[is]

    #p_coes[1,:]=p_coes[1,:][is]
    #p_coes[2,:]=p_coes[2,:][is]
 
    #s_coes[1,:]=s_coes[1,:][is]
    #s_coes[2,:]=s_coes[2,:][is]
    

    ϕ = HelmholtzPotential(bearing.medium.cp, ω / bearing.medium.cp, p_coes, modes)
    ψ = HelmholtzPotential(bearing.medium.cs, ω / bearing.medium.cs, s_coes, modes)

    method = ModalMethod(modes=modes, mode_errors = zeros(T, length(modes)))
    wave = ElasticWave(ω, bearing.medium, [ϕ, ψ], method)
    
    r = boundarytype.inner ? bearing.inner_radius : bearing.outer_radius
   # xs = [r*[cos(θ), sin(θ)] for θ in θs]
    #data =0.0* hcat([field(wave, x, boundarytype.fieldtype) for x in xs]...) |> transpose

    #coefficients = field_modes(wave, r, boundarytype.fieldtype)
    coefficients = field_modes(wave, r, boundarytype.fieldtype, boundarytype)

    return BoundaryData(boundarytype; θs=θs, fields = fouriermodes_to_fields(θs,coefficients,modes), coefficients = coefficients, modes = modes)
end