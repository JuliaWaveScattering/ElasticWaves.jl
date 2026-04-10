
function regular_basis_function_2(medium::Elastic{3,T}, ω::T, field_type::FieldType = DisplacementType()) where T

    return function (order::Integer, x::AbstractVector{T})
        pbasis = pressure_regular_basis(ω, x, medium, order, field_type) 
        Φbasis = shearΦ_regular_basis(ω, x, medium, order, field_type) 
        χbasis = shearχ_regular_basis(ω, x, medium, order, field_type) 

        # the order in each row is first iteration over p, Φ, χ, then increase basis order
        return [pbasis; Φbasis; χbasis] |> transpose
    end
end

function outgoing_basis_function_2(medium::Elastic{3,T}, ω::T, field_type::FieldType = DisplacementType()) where T

    return function (order::Integer, x::AbstractVector{T})
        pbasis = pressure_outgoing_basis(ω, x, medium, order, field_type)
        Φbasis = shearΦ_outgoing_basis(ω, x, medium, order, field_type)
        χbasis = shearχ_outgoing_basis(ω, x, medium, order, field_type)

        return [pbasis; Φbasis; χbasis] |> transpose
    end
end


function t_matrix_2(p::Particle{3,Elastic{3,T},Sphere{T,3}}, outer_medium::Elastic{3,T}, ω::T, basis_order::Integer) where T <: AbstractFloat

    MGφΦs, MGχs = modal_system(p, outer_medium, ω, basis_order) 

    function TmatφΦ(l)
        TφΦ = (MGφΦs[l+1][1] \ MGφΦs[l+1][2])[[1,3],:]
    end

    function Tmatχ(l)
        Tχ  = (MGχs[l+1][1] \ MGχs[l+1][2])[1]
    end

    function Tmatφφ0()
        Tφφ = (MGφΦs[1][1][1:2,1:2] \ MGφΦs[1][2][1:2,1])[1]
        # return [Tφφ zeros(T,1,2); zeros(T,2,3)]
    end

    len(order::Int) = basisorder_to_basislength(PhysicalMedium{3,1},order)

    Tχs = vcat([repeat([Tmatχ(l)],len(l)-len(l-1)) for l in 1:basis_order]...)
    blockcorner = Diagonal(Tχs)

    TmatφΦs = TmatφΦ.(1:basis_order)

    Ts_arr = map(CartesianIndices(TmatφΦs[1])) do inds
        Ts = map(1:basis_order) do l
            repeat([TmatφΦs[l][inds]],len(l)-len(l-1))
        end
        Ts = vcat(Ts...)
    end
    L = length(Ts_arr[1])
    TφΦs_mat = reshape(Ts_arr,2,2)
    DφΦs_mat = Diagonal.(TφΦs_mat)

    Ds = AbstractMatrix{Complex{T}}[Zeros{Complex{T}}(L,L) for i = 1:3, j = 1:3]
    Ds[1:2,1:2] .= DφΦs_mat
    Ds[3,3] = blockcorner

    blocks = map(Ds) do D
        block = Matrix{AbstractMatrix{ComplexF64}}(undef, 2, 2)
        block[1,1] = reshape([zero(Complex{T})],1,1)
        block[1,2] = Zeros{Complex{T}}(1, L)
        block[2,1] = Zeros{Complex{T}}(L,1)
        block[2,2] = D

        return mortar(block)
    end
    # Only the pressure wave has a monopole term, so only the top left corner of the first block is nonzero.
    blocks[1][1,1] = Tmatφφ0()

    return mortar(blocks)
end


"""
    internal_matrix(Particle{3,Elastic{3,T},Sphere{T,3}}, Elastic{3,T}, ω, order)

Similar to the `T-matrix`, except returns the coefficients of the internal field expanded in a basis of regular spherical modes.
"""
function internal_matrix_2(p::Particle{3,Elastic{3,T},Sphere{T,3}}, outer_medium::Elastic{3,T}, ω::T, basis_order::Integer) where T <: AbstractFloat

    MGφΦs, MGχs = modal_system(p, outer_medium, ω, basis_order) 

    function ImatφΦ(l)
        TφΦ = (MGφΦs[l+1][1] \ MGφΦs[l+1][2])[[2,4],:]
    end

    function Imatχ(l)
        Tχ  = (MGχs[l+1][1] \ MGχs[l+1][2])[2]
    end

    function Imatφφ0()
        Tφφ = (MGφΦs[1][1][1:2,1:2] \ MGφΦs[1][2][1:2,1])[2]
    end

    len(order::Int) = basisorder_to_basislength(PhysicalMedium{3,1},order)

    Tχs = vcat([repeat([Imatχ(l)],len(l)-len(l-1)) for l in 1:basis_order]...)
    blockcorner = Diagonal(Tχs)

    TmatφΦs = ImatφΦ.(1:basis_order)

    Ts_arr = map(CartesianIndices(TmatφΦs[1])) do inds
        Ts = map(1:basis_order) do l
            repeat([TmatφΦs[l][inds]],len(l)-len(l-1))
        end
        Ts = vcat(Ts...)
    end
    L = length(Ts_arr[1])
    TφΦs_mat = reshape(Ts_arr,2,2)
    DφΦs_mat = Diagonal.(TφΦs_mat)

    Ds = AbstractMatrix{Complex{T}}[Zeros{Complex{T}}(L,L) for i = 1:3, j = 1:3]
    Ds[1:2,1:2] .= DφΦs_mat
    Ds[3,3] = blockcorner

    blocks = map(Ds) do D
        block = Matrix{AbstractMatrix{ComplexF64}}(undef, 2, 2)
        block[1,1] = reshape([zero(Complex{T})],1,1)
        block[1,2] = Zeros{Complex{T}}(1, L)
        block[2,1] = Zeros{Complex{T}}(L,1)
        block[2,2] = D

        return mortar(block)
    end
    # Only the pressure wave has a monopole term, so only the top left corner of the first block is nonzero.
    blocks[1][1,1] = Imatφφ0()

    return mortar(blocks)
end


ω = 0.9
    basis_order = 10
    field_type = TractionType()

    order = basis_order
    T = Float64

    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)

    centre = [3.0, -3.0, 5.0]
    centre = [0.0, 0.0, 0.0]

    particle_medium = Elastic(3; ρ = 0.6, cp = 2.6, cs = 2.7 ./ 1.2)
    particle_shape = Sphere(centre,1.0)
    particle = Particle(particle_medium, particle_shape)

    order = basis_order;

    # pcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 
    # Φcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 
    # χcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 

    function sourceΦ_coes2(order,centre,ω)
        return [pcoefs Φcoefs χcoefs]
    end
    
    source_field_2 = function (x1,ω) 
        source_basis_2 = regular_basis_function_2(medium, ω, field_type)
        source_basis_2(basis_order, x1 - centre) * sourceΦ_coes2(basis_order,centre,ω)[:] 
    end
    
    sourceΦ_2 = RegularSource{Elastic{3,T},WithoutSymmetry{3}}(medium, source_field_2, sourceΦ_coes2)

    regular_coefficients = regular_spherical_coefficients(sourceΦ_2)
    source_coes_2 = regular_coefficients(basis_order,centre,ω)

    in_matrix = internal_matrix_2(particle, medium, ω, basis_order)
    internal_coes_2 = in_matrix * source_coes_2[:]
    
    scat_matrix = t_matrix_2(particle, medium, ω, basis_order)
    external_coes_2 = scat_matrix * source_coes_2[:]

# choose x on the boundary of the particle
    r = outer_radius(particle)
    xs = [
        centre + spherical_to_cartesian_coordinates([r, i * 2pi / 100, i * 7pi / 100]) 
    for i = 1:100] 
         
    basis = regular_basis_function_2(particle.medium, ω, field_type)
    internal_fields_2 = [basis(basis_order, x - centre) * internal_coes_2[:] for x in xs]
    
    basis = outgoing_basis_function_2(medium, ω, field_type)
    scat_fields_2 = [basis(basis_order, x - centre) * external_coes_2[:] for x in xs]
    source_fields_2 = [source_field_2(x,ω) for x in xs]
    
    external_fields_2 = scat_fields_2 + source_fields_2

    @test norm.(internal_fields_2 - external_fields_2) |> maximum < 1e-13