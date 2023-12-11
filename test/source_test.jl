@testset "Source and scattering" begin

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

