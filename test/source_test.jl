@testset "Source and scattering" begin

    ω = 1.1
    basis_order = 8

    medium = Elastic(3; ρ = 0.5, cp = 1.1, cs = 0.9)

    source = plane_z_shear_source(medium)

    centre = [3.0, -3.0, 5.0]

    source.coefficients(basis_order,centre,ω)
    
    vs = regular_basis_function(source.medium, ω)

    centre = zeros(3)
    
    x = rand(3) .* 0.1 + centre
    
    x = centre

    vs(basis_order, x - centre)
    
    regular_coefficients = regular_spherical_coefficients(source)

    sum(regular_coefficients(basis_order,centre,ω) .* vs(basis_order, x - centre), dims = 1)[:]

    field(source,x,ω) 
    


    source2 = source_expand(source, centre; basis_order = 7)

    xs = [centre + 0.1 .* [cos(τ),sin(τ),rand()] for τ = 0.0:0.3:1.5]
    @test norm([field(source,x,ω) - source2(x,ω) for x in xs]) < 1e-7*norm([field(source,x,ω) for x in xs])


end

