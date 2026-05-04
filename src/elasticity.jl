"""
    Elastic{Dim,T<:AbstractFloat}(ρ::T, c::Complex{T})
    Elastic(ρ::T, c::Union{T,Complex{AbstractFloat}}, Dim::Integer)

Physical properties for a homogenous isotropic elastic medium with wavespeed (c) and density (ρ)

Simulations in this medium produce scalar (Dim) fields in Dim dimensions. In general we use the Debye potentials to describe the field.
"""
struct Elastic{Dim,T} <: PhysicalMedium{Dim,Dim}
    ρ::T # Density (use \greekletter+tab to get the greek letter!)
    cp::Complex{T} # Phase velocity of pressure wave
    cs::Complex{T} # Phase velocity of shear wave
end

# Constructor which supplies the dimension without explicitly mentioning type
function Elastic(Dim::Integer; ρ::T = 0.0, cp::Union{T,Complex{T}} = 0.0, cs::Union{T,Complex{T}} = 0.0) where {T<:Number}
     Elastic{Dim,T}(ρ,Complex{T}(cp),Complex{T}(cs))
end

## wave types

"""
    HelmholtzPotential{T}

See [`field(::HelmholtzPotential, ::AbstractVector)`](@ref) for details on how to evaluate the Helmholtz potential.
"""
struct HelmholtzPotential{T}
    wavespeed::Complex{T}
    wavenumber::Complex{T}
    "The first (second) row is for the besselj (hankelh1) fourier coefficients"
    coefficients::Matrix{Complex{T}}
    modes::Vector{Int}

    function HelmholtzPotential(wavespeed::Complex{T}, wavenumber::Complex{T},  coefficients::AbstractMatrix{Complex{T}}, modes::AbstractVector{Int}) where {T}
        
        # is = sortperm_modes(modes);

        if size(coefficients,1) != 2
            @error "the number of columns in coefficients has to match the basis_order given. There should also be two rows, one for besselj coefficients and another for hankelh1"
        end

        if imag(wavenumber) < 0
            @warn "It is usual to have a wavenumber with a negative imaginary part. Our convention of the Fourier transform implies that this wave is growing exponentially when propagating forward in time."
        end

        new{T}(wavespeed, wavenumber, coefficients, modes |> collect)
    end
end


"""
    ElasticWave{Dim,T}

A type to store an elastic wave expanded in some modal basis. The only basis used in this package is a regular radial expansion, which is needed for scattering and waves within a sphere or cylinder.

For Dim = 4 the field `potentials` are, in order, the φ pressure, Φ shear,  and χ shear, where Φ and χ are the Debye potentials such that the displacement is given by: `u = ∇φ + ∇x∇x(Φ .* r) + ∇x∇x∇x(χ .* r) ./ kS`, where `kS` is the shear wavenumber and `r` is a vector pointing in the radial direction, that is `r = (x,y,z)` in Cartesian coordinates.
"""
struct ElasticWave{Dim,M,T}
    ω::T
    medium::Elastic{Dim,T}
    potentials::Vector{H} where H <:HelmholtzPotential{T}
    method::M

    function ElasticWave(ω::T, medium::Elastic{Dim,T}, potentials::Vector{H}, method::M = ModalMethod();
            mode_errors = zeros(T, length(potentials[1].modes))
        ) where {Dim,T,H <: HelmholtzPotential{T}, M<:SolutionMethod}

        modes_arr = [p.modes for p in potentials]
        lens = length.(modes_arr)

        modes = intersect(modes_arr...)
        len = length(modes)

        if !(findall(lens .!= len) |> isempty)
            @error "The modes for all the potentials is expected to be the same."
        end

        if M == ModalMethod 
            if length(method.mode_errors) != length(potentials[1].modes)
                @warn "The length of mode_errors is expected to be the same as the number of modes."
            end    
        end

        new{Dim,M,T}(ω, medium, potentials, method)
    end
end


import MultipleScattering: name, basisorder_to_basislength, basislength_to_basisorder, regular_basis_function, regular_radial_basis, outgoing_basis_function, outgoing_radial_basis, internal_field
import MultipleScattering: outgoing_translation_matrix, regular_translation_matrix

name(a::Elastic{Dim}) where Dim = "$(Dim)D Elastic"

basisorder_to_basislength(::Type{P}, order::Int) where {T, P<:Elastic{3,T}} = 3 * (order+1)^2
basisorder_to_basislength(::Type{P}, order::Int) where {T, P<:Elastic{2,T}} = 2 * (2*order + 1)

basislength_to_basisorder(::Type{P},len::Int) where {T, P<:Elastic{3,T}} = Int(sqrt(len / 3) - 1)
basislength_to_basisorder(::Type{P},len::Int) where {T, P<:Elastic{2,T}} = Int(T(len / 2 - 1) / T(2.0))

function regular_basis_function(medium::Elastic{3,T}, ω::T, field_type::FieldType = DisplacementType()) where T

    return function (order::Integer, x::AbstractVector{T})
        if norm(x - [0.0,0.0,1.0]) < 1e-10
            throw(ArgumentError("The fields at the north pole of the sphere are not well defined in spherical coordinates."))
            return zeros(Complex{T}, 1, basisorder_to_basislength(Elastic{3,T}, order))
        end
        pbasis = pressure_regular_basis(ω, x, medium, order, field_type) 
        Φbasis = shearΦ_regular_basis(ω, x, medium, order, field_type) 
        χbasis = shearχ_regular_basis(ω, x, medium, order, field_type) 

        # the order in each row is first iteration over p, Φ, χ, then increase basis order
        return reshape([pbasis Φbasis χbasis] |> transpose,3,:)
    end
end

function outgoing_basis_function(medium::Elastic{3,T}, ω::T, field_type::FieldType = DisplacementType()) where T

    return function (order::Integer, x::AbstractVector{T})
        if norm(x - [0.0,0.0,1.0]) < 1e-10
            throw(ArgumentError("The fields at the north pole of the sphere are not well defined in spherical coordinates."))
            return zeros(Complex{T}, 1, basisorder_to_basislength(Elastic{3,T}, order))
        end
        pbasis = pressure_outgoing_basis(ω, x, medium, order, field_type)
        Φbasis = shearΦ_outgoing_basis(ω, x, medium, order, field_type)
        χbasis = shearχ_outgoing_basis(ω, x, medium, order, field_type)

        return reshape([pbasis Φbasis χbasis] |> transpose,3,:)
    end
end

"""
    internal_field(x::AbstractVector, p::Particle{Dim,Elastic{Dim,T}},  source::RegularSource{Elastic{Dim,T}}, ω::T, scattering_coefficients::AbstractVector{Complex{T}})

The internal field for an isotropic elastic particle in an isotropic elastic medium.
"""
function internal_field(x::AbstractVector{T}, p::Particle{Dim,Elastic{Dim,T}}, source::RegularSource{Elastic{Dim,T}}, ω::T, scattering_coefficients::AbstractVector{Complex{T}}, field_type::FieldType = DisplacementType()) where {Dim,T}

    if !(x ∈ p)
        @error "Point $x is not inside the particle with shape $(p.shape)"
    end
    if iszero(p.medium.cp) || isinf(abs(p.medium.cp))
        return zero(Complex{T})
    else
        fs = scattering_coefficients

        order = basislength_to_basisorder(Elastic{Dim,T},length(fs))
        r = outer_radius(p)

        t_mat = t_matrix(p, source.medium, ω, order)
        in_mat = internal_matrix(p, source.medium, ω, order)

        internal_coef0 = [in_mat[1] * (t_mat[1] \ fs[1]), zero(Complex{T}), zero(Complex{T})]
        
        internal_coefs = [internal_coef0; in_mat[4:end,4:end] * (t_mat[4:end,4:end] \ fs[4:end])]

        # # need to seperate the l=0 case
        # L = length(fs)
        # internal_coef0 = in_mat[1,1] * (t_mat[1,1] \ fs[1])

        # # and now remove the l=0 cases
        # inds = [1, (order+1)^2 + 1, 2*(order+1)^2 + 1]
        # in_mat = in_mat[setdiff(1:end, inds), setdiff(1:end, inds)]
        # t_mat = t_mat[setdiff(1:end, inds), setdiff(1:end, inds)]
        # fs = fs[setdiff(1:end, inds)]

        # internal_coefs = zeros(Complex{T}, L)
        # internal_coefs[setdiff(1:end, inds)] .= in_mat * (t_mat \ fs)
        # internal_coefs[1] = internal_coef0

        inner_basis = regular_basis_function(p.medium, ω, field_type)

        return inner_basis(order, x-origin(p)) * internal_coefs
    end
end

import Base.show
function show(io::IO, p::Elastic)
    # Print is the style of the first constructor
    write(io, "Elastic($(p.ρ), $(p.cp),  $(p.cs)) with Dim = $(spatial_dimension(p))")
    return
end

function outgoing_basis_function(medium::Elastic{2}, ω::T, ::PotentialType) where {T<:Number}
    return function (order::Integer, x::AbstractVector{T})
        r, θ  = cartesian_to_radial_coordinates(x)
        kp = ω/medium.cp
        ks = ω/medium.cs
        vcat(
            [hankelh1(m,kp*r)*exp(im*θ*m) for m = -order:order],
            [hankelh1(m,ks*r)*exp(im*θ*m) for m = -order:order]
        ) |> transpose
    end
end

function regular_basis_function(medium::Elastic{2}, ω::T, ::PotentialType) where {T<:Number}
    return function (order::Integer, x::AbstractVector{T})
        r, θ  = cartesian_to_radial_coordinates(x)
        kp = ω/medium.cp
        ks = ω/medium.cs
        vcat(
            [besselj(m,kp*r)*exp(im*θ*m) for m = -order:order],
            [besselj(m,ks*r)*exp(im*θ*m) for m = -order:order]
        ) |> transpose
    end
end

outgoing_translation_matrix(medium::Elastic{3}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where T = outgoing_translation_matrix(medium, in_order, out_order, ω, x, DisplacementType())

function regular_translation_matrix(medium::Elastic{3}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T},::DisplacementType) where {T<:Number}

    # define 2 mediums for pressure and shear potentials 
    pmedium = ScalarMedium{T,3}(medium.cp)
    vp = regular_basis_function(pmedium, ω)(in_order + out_order,x)

    smedium = ScalarMedium{T,3}(medium.cs)
    vs = regular_basis_function(smedium, ω)(in_order + out_order + 1,x)
    
    c(ld,md,l,m,l1) = gaunt_coefficient(T,ld,md,l,m,l1,md-m)

    # auxiliary function to calculate the translation of the shear potentials
        λup(l,m) = sqrt((l-m)*(l+m+1))
        λdown(l,m) = λup(l,-m)
        
        α0(l,m) = sqrt((l+m+1)*(l-m+1)/((2l+1)*(2l+3)))
        β0(l,m) = -sqrt((l-m)*(l+m)/((2l-1)*(2l+1)))
        
        αup(l,m) = -sqrt((l+m+2)*(l+m+1)/((2l+1)*(2l+3)))
        βup(l,m) = -sqrt((l-m)*(l-m-1)/((2l-1)*(2l+1)))

        αdown(l,m) = - αup(l,-m)
        βdown(l,m) = - βup(l,-m)


    # the equivalent of the gaunt coefficients for the shear displacement.
    function C(l,m,dl,dm,l1)
        sumC = zero(Complex{T})
        if abs(m + 1) <= l && abs(dm + 1) <= dl
            sumC += λup(l,m) * λup(l,m) * c(l,m+1,dl,dm+1,l1)
        end
        if abs(m - 1) <= l && abs(dm - 1) <= dl
            sumC += λdown(l,m) * λdown(l,m) * c(l,m-1,dl,dm-1,l1)
        end
        if abs(m) <= l && abs(dm) <= dl && dm != 0 && m != 0
            sumC += 2dm * m * c(l,m,dl,dm,l1)
        end
        sumC = sumC / (2dl * (dl+1))

        return sumC
    end

    function B(ld,md,l,m,l1)
        sumB = zero(Complex{T})
        if l == 0 
            return sumB
        else
            if abs(md - 1) <= ld  
                λnd = λdown(ld,md) 
                if abs(m - 1) <= l + 1
                    sumB +=  -λnd * l * αdown(l,m) * c(ld,md-1,l+1,m-1,l1)
                end    
                if abs(m - 1) <= l - 1
                    sumB +=  λnd * (l+1) * βdown(l,m) * c(ld,md-1,l-1,m-1,l1)
                end
            end 
            if abs(md + 1) <= ld
                λnd = λup(ld,md)
                if abs(m + 1) <= l + 1
                    sumB += -λnd * l * αup(l,m) * c(ld,md+1,l+1,m+1,l1)
                end
                if abs(m + 1) <= l - 1
                    sumB += λnd * (l+1) * βup(l,m) * c(ld,md+1,l-1,m+1,l1)
                end
            end
            if abs(md) <= ld && md != 0
                if abs(m) <= l + 1
                    sumB -= 2md * l * α0(l,-m) * c(ld,md,l+1,m,l1)
                end
                if abs(m) <= l - 1
                    sumB += 2md * (l+1) * β0(l,-m) * c(ld,md,l-1,m,l1)
                end
            end
        end
        sumB = sumB * im / (2l * (l+1))

        return sumB
    end
    
    # the addition translation for each displacement potential 
    ind(order::Int) = basisorder_to_basislength(ScalarMedium{T,3},order)
    U_arr = [
        begin
            i1 = ind(abs(l-dl)-1) + 1;
            i2 = ind(l+dl);

            cs = [
                 (m1 == m - dm && iseven(l1)) ? c(l,m,dl,dm,l1) : 0.0im 
            for l1 = abs(l-dl):(l+dl) for m1 = -l1:l1];

            Cs = [ 
                (m1 == m - dm && iseven(l1) && dl > 0) ? C(l,m,dl,dm,l1) : 0.0im 
            for l1 = abs(l-dl):(l+dl) for m1 = -l1:l1];

            minl1 = min(abs(l-dl+1),abs(l+dl-1));
            maxl1 = l+dl+1;
            Bs = [
                (m1 == m - dm && isodd(l1) && dl > 0) ? B(l,m,dl,dm,l1) : 0.0im 
            for l1 = minl1:maxl1 for m1 = -l1:l1];

            Upp = sum(vp[i1:i2] .* cs)
            Uss = sum(vs[i1:i2] .* Cs)

            i1 = ind(abs(minl1)-1) + 1;
            i2 = ind(maxl1);
            Us  = sum(vs[i1:i2] .* Bs)

            [Upp 0.0im 0.0im; 0.0im Uss Us; 0.0im Us Uss]
        end
    for dl = 0:in_order for dm = -dl:dl for l = 0:out_order for m = -l:l];

    Ublocks = reshape(U_arr, ((out_order+1)^2, (in_order+1)^2))
    U = sparse(mortar(Ublocks))

    return U
end

function outgoing_translation_matrix(medium::Elastic{3}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T},::PotentialType) where {T<:Number}

    # define 2 mediums for pressure and shear potentials 
    pmedium = ScalarMedium{T,3}(medium.cp)
    smedium = ScalarMedium{T,3}(medium.cs)

    U_p = outgoing_translation_matrix(pmedium, in_order, out_order, ω, x)
    U_s = outgoing_translation_matrix(smedium, in_order, out_order, ω, x)

    np,mp = size(U_p) 
    ns,ms = size(U_s)

    # we need to use a block array to multiply with the t_matrix.
    Us = Matrix{AbstractMatrix{Complex{T}}}(undef, 3, 3)
    Us[1,1] = U_p
    Us[1,2] = Zeros{Complex{T}}(np, ms)
    Us[1,3] = Zeros{Complex{T}}(np, ms)
    
    Us[2,2] = U_s
    Us[2,1] = Zeros{Complex{T}}(ns, mp)
    Us[2,3] = Zeros{Complex{T}}(ns, ms)

    Us[3,1] = Zeros{Complex{T}}(ns, mp)
    Us[3,2] = Zeros{Complex{T}}(ns, ms)
    Us[3,3] = U_s

    U = sparse(mortar(Us))

    return U

    #previously:
    # return BlockDiagonal([U_p, U_s, U_s])
end

function regular_translation_matrix(medium::Elastic{3}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}

    # define 2 mediums for pressure and shear potentials 
    pmedium = ScalarMedium{Float64,3}(medium.cp)
    smedium = ScalarMedium{Float64,3}(medium.cs)

    V_p = regular_translation_matrix(pmedium, in_order, out_order, ω, x)
    V_s = regular_translation_matrix(smedium, in_order, out_order, ω, x)

    np,mp = size(V_p)
    ns,ms = size(V_s)

    # we need to use a block array to multiply with the t_matrix.
    Vs = Matrix{AbstractMatrix{Complex{T}}}(undef, 3, 3)
    Vs[1,1] = V_p
    Vs[1,2] = Zeros{Complex{T}}(np, ms)
    Vs[1,3] = Zeros{Complex{T}}(np, ms)
    
    Vs[2,2] = V_s
    Vs[2,1] = Zeros{Complex{T}}(ns, mp)
    Vs[2,3] = Zeros{Complex{T}}(ns, ms)

    Vs[3,1] = Zeros{Complex{T}}(ns, mp)
    Vs[3,2] = Zeros{Complex{T}}(ns, ms)
    Vs[3,3] = V_s

    return sparse(mortar(Vs))
    # return BlockDiagonal([V_p, V_s, V_s])
end

function outgoing_translation_matrix(medium::Elastic{2}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}

    # define 2 mediums for pressure and shear potentials 
    pmedium = ScalarMedium{Float64,2}(medium.cp)
    smedium = ScalarMedium{Float64,2}(medium.cs)

    U_p = outgoing_translation_matrix(pmedium, in_order, out_order, ω, x)
    U_s = outgoing_translation_matrix(smedium, in_order, out_order, ω, x)

    # translation_vec = outgoing_basis_function(medium, ω)(in_order + out_order, x)
    # order = Int(length(translation_vec)/2)
    
    # translation_vec_p = translation_vec[1:order]
    # translation_vec_s = translation_vec[order+1:end]
    # U_p = [
    #     translation_vec_p[n-m + in_order + out_order + 1]
    # for n in -out_order:out_order, m in -in_order:in_order]

    # U_s = [
    #     translation_vec_s[n-m + in_order + out_order + 1]
    # for n in -out_order:out_order, m in -in_order:in_order]

    U = BlockDiagonal([U_p, U_s])
    return U
end

function regular_translation_matrix(medium::Elastic{2}, in_order::Integer, out_order::Integer, ω::T, x::AbstractVector{T}) where {T<:Number}

    # define 2 mediums for pressure and shear potentials 
    pmedium = ScalarMedium{Float64,2}(medium.cp)
    smedium = ScalarMedium{Float64,2}(medium.cs)

    V_p = regular_translation_matrix(pmedium, in_order, out_order, ω, x)
    V_s = regular_translation_matrix(smedium, in_order, out_order, ω, x)

    # translation_vec = regular_basis_function(medium, ω)(in_order + out_order, x)
    # order = Int(length(translation_vec)/2)
    # translation_vec_p = translation_vec[1:order]
    # translation_vec_s = translation_vec[order+1:end]
    # V_p = [
    #     translation_vec_p[n-m + in_order + out_order + 1]
    # for n in -out_order:out_order, m in -in_order:in_order]

    # V_s = [
    #     translation_vec_s[n-m + in_order + out_order + 1]
    # for n in -out_order:out_order, m in -in_order:in_order]

    V = BlockDiagonal([V_p, V_s])
    return V
end