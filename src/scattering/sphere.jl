"""
    t_matrix(Particle{3,Elastic{3,T},Sphere{T,3}}, Elastic{3,T}, ω, order)

The T-matrix for a 3D spherical elastic particle in a 3D elastic medium. This t-matrix `T` has the diagonal matrices of the form `T = diag(T1,T2,...)`. Let `gn = [gφn; gΦn; gχn]` be the coefficients of a regular spherical expansion for the pressure potential φ, and two Debye potentials Φ and χ. See [ElasticSphereScattering.pdf](../../docs/theory/ElasticSphereScattering.pdf) for details. Multiply the matrix with the regular coefficients `Tn * gn` returns the scattering coefficients for the potentials.
"""
function t_matrix(p::Particle{3,Elastic{3,T},Sphere{T,3}}, outer_medium::Elastic{3,T}, ω::T, basis_order::Integer)::Diagonal{Complex{T}} where T <: AbstractFloat

    @warn "There is currently no unit test to show that this T-matrix satisfies the boundary conditions."

    # Sphere radius 
    a = outer_radius(p)

    # wavenumbers of the sphere
    kp1 = ω / p.medium.cp
    ks1 = ω / p.medium.cs

    # wavenumbers of the medium outside the sphere
    kp2 = ω / outer_medium.cp
    ks2 = ω / outer_medium.cs

    # pre-compute spherical bessel and hankel functions
    jkp1 = [sbesselj(l, kp1*a) for l = 0:(basis_order+1)]
    jks1 = [sbesselj(l, ks1*a) for l = 0:(basis_order+1)]
    djkp1 = [diffsbesselj(l, kp1*a) for l = 0:basis_order]
    djks1 = [diffsbesselj(l, ks1*a) for l = 0:basis_order]

    ddjkp1 = [
        2jkp1[l+2] / (kp1 * a) + (l^2 - l - (kp1 *a)^2) * jkp1[l+1] / (kp1 *a)^2
    for l = 0:basis_order]

    hkp1 = [shankelh1(l, kp1*a) for l = 0:(basis_order+1)]
    hks1 = [shankelh1(l, ks1*a) for l = 0:(basis_order+1)]
    dhkp1 = [diffshankelh1(l, kp1*a) for l = 0:basis_order]
    dhks1 = [diffshankelh1(l, ks1*a) for l = 0:basis_order]
    
    jkp2 = [sbesselj(l, kp2*a) for l = 0:(basis_order+1)]
    jks2 = [sbesselj(l, ks2*a) for l = 0:(basis_order+1)]
    djkp2 = [diffsbesselj(l, kp2*a) for l = 0:basis_order]
    djks2 = [diffsbesselj(l, ks2*a) for l = 0:basis_order]

    ddjkp2 = [
        2jkp2[l+2] / (kp2 * a) + (l^2 - l - (kp2 *a)^2) * jkp2[l+1] / (kp2 *a)^2
    for l = 0:basis_order]
    
    hkp2 = [shankelh1(l, kp2*a) for l = 0:(basis_order+1)]
    hks2 = [shankelh1(l, ks2*a) for l = 0:(basis_order+1)]
    dhkp2 = [diffshankelh1(l, kp2*a) for l = 0:basis_order]
    dhks2 = [diffshankelh1(l, ks2*a) for l = 0:basis_order]

    ddjs = [
        2js[l+2] / ksr + (l^2 - l - ksr^2) * js[l+1] / ksr^2
    for l = 0:basis_order]

    function Mχ(l::Integer)
        a * ks2 

    "Returns a ratio used in multiple scattering which reflects the material properties of the particles"
    function Zn(m::Integer)::Complex{T}
        m = T(abs(m))
        ak = outer_radius(p)*ω/outer_medium.c

        # set the scattering strength and type
        if isinf(p.medium.c) || isinf(p.medium.ρ) || iszero(outer_medium.ρ)
            numer = diffsbesselj(m, ak)
            denom = diffshankelh1(m, ak)
        elseif iszero(p.medium.ρ) || iszero(p.medium.c)
            numer = sbesselj(m, ak)
            denom = shankelh1(m, ak)
        else
            q = impedance(p.medium)/impedance(outer_medium) # Impedance ratio
            γ = outer_medium.c/p.medium.c #speed ratio
            numer = q * diffsbesselj(m, ak) * sbesselj(m, γ * ak) - sbesselj(m, ak) * diffsbesselj(m, γ * ak)
            denom = q * diffshankelh1(m, ak) * sbesselj(m, γ * ak) - shankelh1(m, ak) * diffsbesselj(m, γ * ak)
        end

        return numer / denom
    end
    Zns = Zn.(0:basis_order)

    len(order::Int) = basisorder_to_basislength(Acoustic{T,3},order)
    T_vec = - vcat(Zns[1],[repeat(Zns[l+1:l+1],len(l)-len(l-1)) for l = 1:basis_order]...)

    return Diagonal(T_vec)
end
