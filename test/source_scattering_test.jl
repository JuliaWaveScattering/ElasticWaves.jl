@testset "Source" begin

    ω = 1.2
    basis_order = 16
    field_type = DisplacementType()
    order = basis_order

    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)
    ks = ω / medium.cs

    source = plane_z_shear_source(medium)

    centre = [3.0, -3.0, 5.0]
    x = rand(3) .* 0.5 + centre
    
    regular_coefficients = regular_spherical_coefficients(source)
    coes = regular_coefficients(basis_order,centre,ω)
    
    basis = regular_basis_function(medium, ω)
    field_reg = basis(basis_order, x - centre) * regular_coefficients(basis_order,centre,ω)[:] 

    @test field_reg - field(source,x,ω) |> norm < 2e-14


end

@testset "Source" begin

    ω = 1.2
    basis_order = 16
    field_type = DisplacementType()
    order = basis_order

    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)
    ks = ω / medium.cs

    source = plane_z_shear_source(medium)

    centre = [3.0, -3.0, 5.0]
    x = rand(3) .* 0.5 + centre
    
    regular_coefficients = regular_spherical_coefficients(source)
    coes = regular_coefficients(basis_order,centre,ω)
    
    basis = regular_basis_function(medium, ω)
    field_reg = basis(basis_order, x - centre) * regular_coefficients(basis_order,centre,ω)[:] 

    @test field_reg - field(source,x,ω) |> norm < 2e-14

    particle_medium = Elastic(3; ρ = 0.6, cp = 0.6, cs = 0.7 ./ 1.2)
    particle_shape = Sphere(3,1.0)
    particle = Particle(particle_medium, particle_shape)

    # use MultipleScattering to calculate the scattered wave.
    t_matrix(particle, medium, ω, basis_order)


end

