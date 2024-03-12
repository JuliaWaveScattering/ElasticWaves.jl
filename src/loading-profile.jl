# Here are all functions and possibly types related to the loading of the roller bearing raceway 

function natural_frequencies(bearing::RollerBearing, frequency_order::Int)
    Ω = bearing.angular_speed
    Z = bearing.number_of_rollers

    ωs = (0:frequency_order) .* (Z * Ω)

    return ωs
end    


"""
    BoundaryData(ω::AbstractFloat, bearing::RollerBearing, loading_profile::BoundaryData)

    Return BoundaryData for bearing when assume point contact of the roller bearings and using the loading_profile (what load the bearings hold).

    The frequency of the loading profile data provided is assumed to be ω - ω_m, where m = round(ω / (Z * Ω)). This is a bit opaque, so in the future we need to rewrite this to make it clearer.
"""    
function BoundaryData(ω::Number, bearing::RollerBearing, loading_profile::BoundaryData)

    Ω = bearing.angular_speed
    Z = bearing.number_of_rollers

    if isempty(loading_profile.fourier_modes)
        basis_order = Int(floor(length(loading_profile.θs)/2.0 - 1/2.0))
        loading_profile = fields_to_fouriermodes(loading_profile, basis_order)
    end
        
    fourier_modes = loading_profile.fourier_modes 
    
    basis_order = basislength_to_basisorder(PhysicalMedium{2,1}, size(fourier_modes,1))

    # the natural wavenumber of the bearing ω_m that is closest to the given frequency ω is ω_m = frequency_order * Z * Ω
    frequency_order = Int(round(ω / (Z * Ω)))

    fourier_order = 2*(basis_order + Z*frequency_order) + 1 

    bd_fourier_modes = zeros(ComplexF64,fourier_order,2) 
    bd_fourier_modes[2*Z*frequency_order+1:end,:] = (Z/(2pi)) .* fourier_modes

    fields = fouriermodes_to_fields(loading_profile.θs, fourier_modes)

    bd =  BoundaryData(loading_profile.boundarytype; 
        θs = loading_profile.θs, 
        fields = fields,
        fourier_modes = bd_fourier_modes
    )

    return bd
end

function point_contact_boundary_data(θs::Vector{T}, bearing::RollerBearing, bc ::BoundaryCondition{TractionType}, basis_order::Int, friction_coefficient = 1.0) where{T}

    Z =bearing.number_of_rollers
    # μ=bearing.medium.friction_coefficient
    n_order=basis_order

    size=2*n_order*Z+1 

    fourier_coef_p=zeros(size)

    for n in -n_order:n_order
         fourier_coef_p[Int(n_order*Z)+1+Int(n*Z)]=Z/2pi
    end

    fourier_coef_s=μ.*fourier_coef_p 
    bd =  BoundaryData(bc, θs=θs, fourier_modes=hcat(fourier_coef_p,fourier_coef_s))

    return bd
    
end
