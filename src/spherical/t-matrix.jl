import MultipleScattering: t_matrix


"""
    modal_system(Particle{3,Elastic{3,T},Sphere{T,3}}, Elastic{3,T}, ω, order)

Returns an array of the modal matrices a 3D spherical elastic particle in a 3D elastic medium. See `t_matrix` for more details.
"""
function modal_system(p::Particle{3,Elastic{3,T},Sphere{T,3}}, outer_medium::Elastic{3,T}, ω::T, basis_order::Integer) where T <: AbstractFloat

    # Sphere radius 
    a = outer_radius(p)

    # wavenumbers of the sphere
    kp1 = ω / p.medium.cp
    ap1 = a*kp1

    ks1 = ω / p.medium.cs
    as1 = a*ks1
    
    # wavenumbers of the medium outside the sphere
    kp2 = ω / outer_medium.cp
    ap2 = a*kp2
    
    ks2 = ω / outer_medium.cs
    as2 = a*ks2

    # the Lame parameters: cs = sqrt(μ/ρ), cp = sqrt((λ + 2μ)/ρ)
    μ1 = p.medium.ρ * p.medium.cs^2
    λ1 = p.medium.cp^2 * p.medium.ρ - 2μ1 
    
    μ2 = outer_medium.ρ * outer_medium.cs^2
    λ2 = outer_medium.cp^2 * outer_medium.ρ - 2μ2 

    # pre-compute spherical bessel and hankel functions
    jp1   = [     sbesselj(l, ap1) for l = 0:basis_order]
    djp1  = [ diffsbesselj(l, ap1) for l = 0:basis_order]
    ddjp1 = [diff2sbesselj(l, ap1) for l = 0:basis_order]
    
    js1    = [     sbesselj(l, as1) for l = 0:basis_order]
    djs1   = [ diffsbesselj(l, as1) for l = 0:basis_order]
    ddjs1  = [diff2sbesselj(l, as1) for l = 0:basis_order]
    dddjs1 = [diff3sbesselj(l, as1) for l = 0:basis_order]
    
    js2    = [     sbesselj(l, as2) for l = 0:basis_order]
    djs2   = [ diffsbesselj(l, as2) for l = 0:basis_order]
    ddjs2  = [diff2sbesselj(l, as2) for l = 0:basis_order]
    dddjs2 = [diff3sbesselj(l, as2) for l = 0:basis_order]

    jp2   = [     sbesselj(l, ap2) for l = 0:basis_order]
    djp2  = [ diffsbesselj(l, ap2) for l = 0:basis_order]
    ddjp2 = [diff2sbesselj(l, ap2) for l = 0:basis_order]

    hs2    = [     shankelh1(l, as2) for l = 0:basis_order]
    dhs2   = [ diffshankelh1(l, as2) for l = 0:basis_order]
    ddhs2  = [diff2shankelh1(l, as2) for l = 0:basis_order]
    dddhs2 = [diff3shankelh1(l, as2) for l = 0:basis_order]

    hp2   = [     shankelh1(l, ap2) for l = 0:basis_order]
    dhp2  = [ diffshankelh1(l, ap2) for l = 0:basis_order]
    ddhp2 = [diff2shankelh1(l, ap2) for l = 0:basis_order]
    
    Mχ(l) = [as2*hs2[l+1]                    -as1*js1[l+1]; 
            as2*μ2*(as2*dhs2[l+1]-hs2[l+1]) as1*μ1*(js1[l+1]-as1*djs1[l+1])] 
        
    Gχ(l) = [-as2*js2[l+1], as2*μ2*(js2[l+1]-as2*djs2[l+1])]      

    Gφ(l) = [ap2*djp2[l+1],
        ap2^2*(2μ2*ddjp2[l+1]-λ2*jp2[l+1]),
        jp2[l+1], 
        2μ2*(ap2*djp2[l+1]-jp2[l+1])
    ]
     
    GΦ(l) = [
        as2*(2*djs2[l+1] + as2*js2[l+1] + as2*ddjs2[l+1]), 
        2as2^2*μ2*(js2[l+1] + as2*djs2[l+1] + 3*ddjs2[l+1] + as2*dddjs2[l+1]), 
        js2[l+1] + as2*djs2[l+1],
        μ2*((as2^2-2)*js2[l+1] + 2*as2*(djs2[l+1] + as2*ddjs2[l+1]))
    ]

    M11(l) = [
        -ap2*dhp2[l+1]                     ap1*djp1[l+1];
        ap2^2*(λ2*hp2[l+1]-2μ2*ddhp2[l+1]) ap1^2*(2μ1*ddjp1[l+1]-λ1*jp1[l+1]) 
    ]

    M21(l) = [
        -hp2[l+1]                     jp1[l+1];
        (2μ2*hp2[l+1]-2μ2*ap2*dhp2[l+1]) (2μ1*ap1*djp1[l+1]-2*μ1*jp1[l+1]) 
    ]

    M12(l) = [
        -as2*(as2*hs2[l+1]+2dhs2[l+1]+as2*ddhs2[l+1])                         as1*(as1*js1[l+1] + 2*djs1[l+1] + as1*ddjs1[l+1]); 
        -2as2^2*μ2*(hs2[l+1] + as2*dhs2[l+1] + 3ddhs2[l+1] + as2*dddhs2[l+1]) 2as1^2*μ1*(js1[l+1] + as1*djs1[l+1] + 3ddjs1[l+1] + as1*dddjs1[l+1]) 
    ]

    M22(l) = [
        (-hs2[l+1]-as2*dhs2[l+1])                                 js1[l+1]+as1*djs1[l+1];
        (μ2*(2-as2^2)*hs2[l+1]-2μ2*as2*(dhs2[l+1]+as2*ddhs2[l+1])) (μ1*(as1^2-2)*js1[l+1]+2μ1*as1*(djs1[l+1]+as1*ddjs1[l+1]))
    ]

    MGφΦs = [
        [[M11(l) M12(l); 
         M21(l) M22(l)],
         [Gφ(l) GΦ(l)]]
    for l = 0:basis_order]
    
    MGχs = [[Mχ(l),Gχ(l)] for l = 0:basis_order]

    return MGφΦs, MGχs
end    

"""
    t_matrix(Particle{3,Elastic{3,T},Sphere{T,3}}, Elastic{3,T}, ω, order)

The T-matrix for a 3D spherical elastic particle in a 3D elastic medium. This t-matrix `T` has the diagonal matrices of the form `T = diag(T1,T2,...)`. Let `gn = [gφn; gΦn; gχn]` be the coefficients of a regular spherical expansion for the pressure potential φ, and two Debye potentials Φ and χ. See [ElasticSphereScattering.pdf](../../docs/theory/ElasticSphereScattering.pdf) for details. Multiply the matrix with the regular coefficients `Tn * gn` returns the scattering coefficients for the potentials.
"""
function t_matrix(p::Particle{3,Elastic{3,T},Sphere{T,3}}, outer_medium::Elastic{3,T}, ω::T, basis_order::Integer) where T <: AbstractFloat

    # @warn "There is a unit test to show that this T-matrix satisfies displacement boundary conditions but have not yet tested the traction boundary."
    
    MGφΦs, MGχs = modal_system(p, outer_medium, ω, basis_order) 

    function Tmat(l)
        TφΦ = (MGφΦs[l+1][1] \ MGφΦs[l+1][2])[[1,3],:]
        Tχ  = (MGχs[l+1][1] \ MGχs[l+1][2])[1]
        return [TφΦ zeros(T,2); T(0) T(0) Tχ]
    end

    function Tmat0()
        Tφφ = (MGφΦs[1][1][1:2,1:2] \ MGφΦs[1][2][1:2,1])[1]
        return [Tφφ zeros(T,1,2); zeros(T,2,3)]
    end

    Tmats = Tmat.(1:basis_order)

    len(order::Int) = basisorder_to_basislength(PhysicalMedium{3,1},order)
    T_vec = [Tmat0(),
        vcat([repeat(Tmats[l:l],len(l)-len(l-1)) for l = 1:basis_order]...)...
    ]

    return BlockDiagonal(T_vec)
end

"""
    internal_matrix(Particle{3,Elastic{3,T},Sphere{T,3}}, Elastic{3,T}, ω, order)

Similar to the `T-matrix`, except returns the coefficients of the internal field expanded in a basis of regular spherical modes.
"""
function internal_matrix(p::Particle{3,Elastic{3,T},Sphere{T,3}}, outer_medium::Elastic{3,T}, ω::T, basis_order::Integer) where T <: AbstractFloat

    MGφΦs, MGχs = modal_system(p, outer_medium, ω, basis_order) 

    function inner_mat(l)
        TφΦ = (MGφΦs[l+1][1] \ MGφΦs[l+1][2])[[2,4],:]
        Tχ  = (MGχs[l+1][1] \ MGχs[l+1][2])[2]
        return [TφΦ zeros(T,2); T(0) T(0) Tχ]
    end

    function inner_mat0()
        Tφφ = (MGφΦs[1][1][1:2,1:2] \ MGφΦs[1][2][1:2,1])[2]
        return [Tφφ zeros(T,1,2); zeros(T,2,3)]
    end

    inner_mats = inner_mat.(1:basis_order)

    len(order::Int) = basisorder_to_basislength(PhysicalMedium{3,1},order)
    T_vec = [inner_mat0(),
        vcat(
            [repeat(inner_mats[l:l],len(l)-len(l-1)) for l = 1:basis_order]...
        )...
    ]

    return BlockDiagonal(T_vec)
end
