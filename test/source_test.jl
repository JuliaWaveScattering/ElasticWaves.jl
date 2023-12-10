@testset "Source and scattering" begin

    ω = 1.1
    basis_order = 12
    field_type = DisplacementType()
    order = basis_order

    medium = Elastic(3; ρ = 0.5, cp = 1.1, cs = 0.9)

    source = plane_z_shear_source(medium)

    centre = [3.0, -3.0, 5.0]


    basis = regular_basis_function(medium, ω)

    centre = zeros(3)
    
    x = rand(3) .* 0.1 + centre
    x = rand(3) .* 0.000001 + centre
    
    # basis(basis_order, x - centre)
    regular_coefficients = regular_spherical_coefficients(source)

    field_reg = basis(basis_order, x - centre) * regular_coefficients(basis_order,centre,ω)[:] 

    rθφ = cartesian_to_spherical_coordinates(x - centre)
    
    M1 = spherical_to_cartesian_transform(rθφ)  
    M1 * field_reg
    
    M2 = cartesian_to_spherical_transform(x - centre)
    M2 * M1
    field(source,x,ω) 
    


    # source2 = source_expand(source, centre; basis_order = 7)

    # xs = [centre + 0.1 .* [cos(τ),sin(τ),rand()] for τ = 0.0:0.3:1.5]
    # @test norm([field(source,x,ω) - source2(x,ω) for x in xs]) < 1e-7*norm([field(source,x,ω) for x in xs])


end

