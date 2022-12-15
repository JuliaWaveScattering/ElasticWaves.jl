function boundarycondition_modes(ω::AbstractFloat, basis_order::Int, bc::TractionBoundary, bearing::Bearing{2})
    r = (bc.inner == true) ? bearing.r1 : bearing.r2
    return [
        pressure_traction_modes(ω, r, basis_order, bearing.medium);
        shear_traction_modes(ω, r, basis_order, bearing.medium)
    ]
end

function boundarycondition_modes(ω::AbstractFloat, basis_order::Int, bc::DisplacementBoundary, bearing::Bearing{2})
    r = (bc.inner == true) ? bearing.r1 : bearing.r2
    return [
        pressure_displacement_modes(ω, r, basis_order, bearing.medium);
        shear_displacement_modes(ω, r, basis_order, bearing.medium)
    ]
end


function boundarycondition_mode(ω::AbstractFloat, basis_order::Int, bc::TractionBoundary, bearing::Bearing{2})
    r = (bc.inner == true) ? bearing.r1 : bearing.r2
    return hcat(pressure_traction_mode(ω, r, basis_order, bearing.medium), shear_traction_mode(ω, r, basis_order, bearing.medium))
end

function boundarycondition_mode(ω::AbstractFloat, basis_order::Int, bc::DisplacementBoundary, bearing::Bearing{2})
    r = (bc.inner == true) ? bearing.r1 : bearing.r2
    return hcat(pressure_displacement_mode(ω, r, basis_order, bearing.medium), shear_displacement_mode(ω, r, basis_order, bearing.medium))
end

function boundarycondition_system(ω::AbstractFloat, basis_order::Int, bearing::Bearing{2}, bcs::AbstractBoundaryConditions)

    first_boundary = boundarycondition_mode(ω, basis_order, bcs[1], bearing)
    second_boundary = boundarycondition_mode(ω, basis_order, bcs[2], bearing)

    return vcat(first_boundary, second_boundary)
end

function boundarycondition_systems(ω::AbstractFloat, basis_order::Int, bearing::Bearing{2}, bcs::AbstractBoundaryConditions)

    p1_j, p1_h, s1_j, s1_h = boundarycondition_modes(ω, basis_order, bcs[1], bearing);
    p2_j, p2_h, s2_j, s2_h = boundarycondition_modes(ω, basis_order, bcs[2], bearing);

    return M = map(1:(2basis_order+1)) do n
        vcat(
            hcat(p1_j[:,n],p1_h[:,n],s1_j[:,n],s1_h[:,n]),
            hcat(p2_j[:,n],p2_h[:,n],s2_j[:,n],s2_h[:,n])
        );
    end
end

function ElasticWave(ω::T, bearing::Bearing{2}, bcs::AbstractBoundaryConditions, forcing::Matrix{C}) where {T, C <: Complex{T}}

    kP = ω / bearing.medium.cp;
    kS = ω / bearing.medium.cs

    basis_order = Int( (size(forcing,1) - 1) / 2)

    coefficients = map(-basis_order:basis_order) do m
        A = boundarycondition_system(ω, m, bearing, bcs)
        b = forcing[m+basis_order+1,:]
        A \ b
    end

    coefficients = transpose(hcat(coefficients...)) |> collect

    pressure_bessel_coefficients = coefficients[:,1]
    pressure_hankel_coefficients = coefficients[:,2]

    shear_bessel_coefficients = coefficients[:,3]
    shear_hankel_coefficients = coefficients[:,4]

    φ = HelmholtzPotential{2}(kP,basis_order,pressure_bessel_coefficients,pressure_hankel_coefficients)
    ψ = HelmholtzPotential{2}(kS,basis_order,shear_bessel_coefficients,shear_hankel_coefficients)

    return ElasticWave{2}(ω,bearing.medium, φ, ψ)
end
