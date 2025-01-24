# Create elastic wave sources
import MultipleScattering: point_source, AbstractSource, outgoing_basis_function, regular_basis_function, outgoing_translation_matrix, regular_translation_matrix

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
        # p_potential = HelmholtzPotential{3}(medium.cp, ω / medium.cp, [pcoefs; 0 .* pcoefs])
        
        Φcoefs = T(sqrt(pi)) * sum(source_field(centre,ω) .* polarisation) .*
        [
            (abs(m) == 1) ?  Complex{T}(m * (1.0im)^l * sqrt(T(2l + 1) / T((1+l)*l))) : Complex{T}(0)
        for l = 0:order for m = -l:l] 
        # Φ_potential = HelmholtzPotential{3}(medium.cs, ω / medium.cs, [Φcoefs; 0 .* Φcoefs])
        
        χcoefs = T(sqrt(pi)) * sum(source_field(centre,ω) .* polarisation) .*
        [
            (abs(m) == 1) ?  Complex{T}((1.0im)^l * sqrt(T(2l + 1) / T((1+l)*l))) : Complex{T}(0)
        for l = 0:order for m = -l:l] 
        # χ_potential = HelmholtzPotential{3}(medium.cs, ω / medium.cs, [χcoefs; 0 .* χcoefs])

        # return ElasticWave(ω, medium, [pcoefs,Φcoefs,χcoefs])
        return [pcoefs Φcoefs χcoefs] ./ (-im * ks) |> transpose
    end

    return RegularSource{Elastic{3,T},S}(medium, source_field, spherical_expansion)
end

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

function outgoing_basis_function(medium::Elastic{2}, ω::T) where {T<:Number}
    return function (order::Integer, x::AbstractVector{T})
        r, θ  = cartesian_to_radial_coordinates(x)
        kp = ω/medium.cp
        ks = ω/medium.cs
        vcat([hankelh1(m,kp*r)*exp(im*θ*m) for m = -order:order],[hankelh1(m,ks*r)*exp(im*θ*m) for m = -order:order]) |> transpose
    end
end

function regular_basis_function(medium::Elastic{2}, ω::T) where {T<:Number}
    return function (order::Integer, x::AbstractVector{T})
        r, θ  = cartesian_to_radial_coordinates(x)
        kp = ω/medium.cp
        ks = ω/medium.cs
        vcat([besselj(m,kp*r)*exp(im*θ*m) for m = -order:order],[besselj(m,ks*r)*exp(im*θ*m) for m = -order:order]) |> transpose
    end
end

function outgoing_translation_matrix(medium::Elastic{2}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}
    translation_vec = outgoing_basis_function(medium, ω)(in_order + out_order, x)
    order = Int(length(translation_vec)/2)
    translation_vec_p = translation_vec[1:order]
    translation_vec_s = translation_vec[order+1:end]
    U_p = [
        translation_vec_p[n-m + in_order + out_order + 1]
    for n in -out_order:out_order, m in -in_order:in_order]

    U_s = [
        translation_vec_s[n-m + in_order + out_order + 1]
    for n in -out_order:out_order, m in -in_order:in_order]

    U = BlockDiagonal([U_p, U_s])
    return U
end

function regular_translation_matrix(medium::Elastic{2}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}
    translation_vec = regular_basis_function(medium, ω)(in_order + out_order, x)
    order = Int(length(translation_vec)/2)
    translation_vec_p = translation_vec[1:order]
    translation_vec_s = translation_vec[order+1:end]
    V_p = [
        translation_vec_p[n-m + in_order + out_order + 1]
    for n in -out_order:out_order, m in -in_order:in_order]

    V_s = [
        translation_vec_s[n-m + in_order + out_order + 1]
    for n in -out_order:out_order, m in -in_order:in_order]

    V = BlockDiagonal([V_p, V_s])
    return V
end


function boundary_data(ω::Float64, bearing::RollerBearing, boundarytype::BC, p_map::SourceMap, s_map::SourceMap, modes::Vector{Int}, θs::AbstractVector) where {BC <: BoundaryCondition}
    fieldtype = boundarytype.fieldtype
    inner = boundarytype.inner
    order = maximum(modes)
    basis_length = 2*order+1

    locations_p = p_map.locations
    amplitudes_p = p_map.amplitudes
    locations_s = s_map.locations
    amplitudes_s = s_map.amplitudes

    number_of_p_sources = length(amplitudes_p)
    number_of_s_sources = length(amplitudes_s)

    r = inner ? bearing.inner_radius : bearing.outer_radius

    xs = [r*[cos(θ), sin(θ)] for θ in θs]
    
    if inner == true
        # potentials of the form ΣH(kr')J(kr)exp(i*n*(θ-θ'))

        p_coes = [im*amplitudes_p[i]*hcat(outgoing_translation_matrix(bearing.medium, order, 0, ω, -locations_p[i,:])[1,:][1:basis_length], zeros(ComplexF64, basis_length))/4 for i in 1:number_of_p_sources]
        p_coes = sum(p_coes)
        s_coes = [im*amplitudes_s[i]*hcat(outgoing_translation_matrix(bearing.medium, order, 0, ω, -locations_s[i,:])[2,:][basis_length+1:end], zeros(ComplexF64, basis_length))/4 for i in 1:number_of_s_sources]
        s_coes = sum(s_coes)

        ϕ = HelmholtzPotential{2}(bearing.medium.cp, ω/bearing.medium.cp, p_coes |> transpose, modes)
        ψ = HelmholtzPotential{2}(bearing.medium.cs, ω/bearing.medium.cs, s_coes |> transpose, modes)

        wave = ElasticWave(ω, bearing.medium, [ϕ, ψ])
        
        data = hcat([field(wave, x, fieldtype) for x in xs]...) |> transpose
    else
        # potentials of the form ΣJ(kr')H(kr)exp(i*n*(θ-θ'))

        p_coes = [im*amplitudes_p[i]*hcat(zeros(ComplexF64, basis_length),regular_translation_matrix(bearing.medium, order, 0, ω, -locations_p[i,:])[1,:][1:basis_length])/4 for i in 1:number_of_p_sources]
        p_coes = sum(p_coes)
        
        s_coes = [im*amplitudes_s[i]*hcat(zeros(ComplexF64, basis_length),regular_translation_matrix(bearing.medium, order, 0, ω, -locations_s[i,:])[2,:][basis_length+1:end])/4 for i in 1:number_of_s_sources]
        s_coes = sum(s_coes)

        ϕ = HelmholtzPotential{2}(bearing.medium.cp, ω/bearing.medium.cp, p_coes |> transpose, modes)
        ψ = HelmholtzPotential{2}(bearing.medium.cs, ω/bearing.medium.cs, s_coes |> transpose, modes)

        wave = ElasticWave(ω, bearing.medium, [ϕ, ψ])
        
        data = hcat([field(wave, x, fieldtype) for x in xs]...) |> transpose
    end

    return BoundaryData(boundarytype; θs=θs, fields = Matrix(data))
end