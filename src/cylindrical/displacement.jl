
function potentialϕ(order::Int;ω::T,r::T,θ::T,bearing::Bearing{T,2},fp_coefficients::Vector,fs_coefficients::Vector) where {T<:Number}
    n = order
    f = coes(n;ω,bearing,fp_coefficients, fs_coefficients)
    kp = ω/bearing.medium.cp
    return (f[1]*besselj(n,kp*r)+f[2]*hankelh1(n,kp*r))*ℯ^(im*n*θ)
end

function potentialψ(order::Int;ω::T,r::T,θ::T,bearing::Bearing{T,2},fp_coefficients::Vector,fs_coefficients::Vector) where {T<:Number}
    n = order
    f = coes(n;ω,bearing,fp_coefficients, fs_coefficients)
    ks = ω/bearing.medium.cs
    return (f[3]*besselj(n,ks*r)+f[4]*hankelh1(n,ks*r))*ℯ^(im*n*θ)
end

function displacement_mode(order::Int;ω::T,r::T,θ::T,bearing::Bearing{T,2},fp_coefficients::Vector,fs_coefficients::Vector) where {T<:Number}
    n = order
    kp = ω/bearing.medium.cp
    ks = ω/bearing.medium.cs
    f = coes(n;ω,bearing,fp_coefficients,fs_coefficients)

    u1 = cos(θ) * (
            kp * (f[1]*diffbesselj(n,kp*r) + f[2]*diffhankelh1(n,kp*r)) +
            ((im*n)/r) * potentialψ(n;ω,r,θ,bearing,fp_coefficients,fs_coefficients)
        ) - sin(θ) * (
            ((im*n)/r) * potentialϕ(n;ω,r,θ,bearing,fp_coefficients,fs_coefficients) -ks*(f[3]*diffbesselj(n,ks*r)+f[4]*diffhankelh1(n,ks*r))
        )

    u2 = sin(θ) * (
            kp * (f[1]*diffbesselj(n,kp*r)+f[2]*diffhankelh1(n,kp*r)) +
            ((im*n)/r) * potentialψ(n;ω,r,θ,bearing,fp_coefficients,fs_coefficients)
        ) + cos(θ) * (
            ((im*n)/r) * potentialϕ(n;ω,r,θ,bearing,fp_coefficients,fs_coefficients) -ks*(f[3]*diffbesselj(n,ks*r)+f[4]*diffhankelh1(n,ks*r))
        )

    return [u1;u2]*ℯ^(im*n*θ)
end

function displacement(basis_order::Int;ω::T,r::T,θ::T,bearing::Bearing{T,2},fp_coefficients::Vector,fs_coefficients::Vector) where {T<:Number}

    u=zeros(2,1)

    for i in 1:basis_order
        u=u+displacement_mode(convert(Int,i-(basis_order+1)/2);ω,r,θ,bearing,fp_coefficients, fs_coefficients)
    end

    return u
end
