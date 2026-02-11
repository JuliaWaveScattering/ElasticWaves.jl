import MultipleScattering: t_matrix

function t_matrix(p::Particle{2,Elastic{2,T},Sphere{T,2}}, outer_medium::Elastic{2,T}, ω::T, basis_order::Integer) where T <: AbstractFloat

    # Circle radius 
    a = outer_radius(p)
    ρ = outer_medium.ρ
    kP=ω / outer_medium.cp
    kS=ω / outer_medium.cs
    modes=-basis_order:basis_order
    cs=outer_medium.cs

    ρs = ρ * cs^2 / a^2
    bessel_modesP(J::Function,n) = [
        -ρs * (2kP * a * J(n-1, kP*a) + (- 2n - 2n^2 + kS^2 * a^2) * J(n, kP*a)),
        im * n * ρs * (kP*a*J(n-1, kP*a) - 2J(n, kP*a) - kP * a * J(1 + n, kP*a))
    ]

    ρs2 = ρ * cs^2 / a^2
    bessel_modesS(J::Function,n) = [
        im * n * ρs2 * (kS*a*J(n-1, kS*a) - 2J(n, kS*a) - kS*a*J(n+1, kS*a)),
        ρs2 * (2kS*a*J(n-1, kS*a) + (-2*n - 2*n^2 + kS^2 * a^2) * J(n, kS*a))
    ]
    
    m13=[hcat(bessel_modesP(besselj,n), bessel_modesS(besselj,n))  for n in modes]
    m24=[hcat(bessel_modesP(hankelh1,n), bessel_modesS(hankelh1,n))  for n in modes]
    Ts=[-inv(m24[i])*m13[i] for i in eachindex(modes)]

    return BlockDiagonal(Ts)

end

function modal_system(p::Particle{2,Elastic{2,T},Sphere{T,2}}, outer_medium::Elastic{2,T}, ω::T, basis_order::Integer) where T <: AbstractFloat

    
    # Cylinder radius 
    a = outer_radius(p)
    ρ = outer_medium.ρ
    kP=ω / outer_medium.cp
    kS=ω / outer_medium.cs
    modes=-basis_order:basis_order
    cs=outer_medium.cs

    ρs = ρ * cs^2 / a^2
    bessel_modesP(J::Function,n) = [
        -ρs * (2kP * a * J(n-1, kP*a) + (- 2n - 2n^2 + kS^2 * a^2) * J(n, kP*a)),
        im * n * ρs * (kP*a*J(n-1, kP*a) - 2J(n, kP*a) - kP * a * J(1 + n, kP*a))
    ]

    ρs2 = ρ * cs^2 / a^2
    bessel_modesS(J::Function,n) = [
        im * n * ρs2 * (kS*a*J(n-1, kS*a) - 2J(n, kS*a) - kS*a*J(n+1, kS*a)),
        ρs2 * (2kS*a*J(n-1, kS*a) + (-2*n - 2*n^2 + kS^2 * a^2) * J(n, kS*a))
    ]
    
    m13=[hcat(bessel_modesP(besselj,n), bessel_modesS(besselj,n))  for n in modes]
    m24=[hcat(bessel_modesP(hankelh1,n), bessel_modesS(hankelh1,n))  for n in modes]
   

    return BlockDiagonal(m13), BlockDiagonal(m24)

    
end    

