@testset "Translation 2D Matrices test" begin

    ω = rand() + 0.1
    medium = Elastic(2; ρ = 1.0, cp = 1.0, cs = 0.5)
    r = rand(2) - [0.5,0.5];
    d = rand(2) - [0.5,0.5];
    d = 10 * d * norm(r) / norm(d)

    # Note that to be accurate the order of vs
    order = 4
    larger_order = 3*order
    basis_length = 2*larger_order+1
    
    # Test 2D outgoing translation matrix
    U = outgoing_translation_matrix(medium, larger_order, order, ω, d)

    field_type = PotentialType()
    vs = regular_basis_function(medium, ω, field_type)(larger_order,r)
    us = outgoing_basis_function(medium, ω, field_type)(order,r + d)

    @test maximum(abs.(U * vs[:] - us[:]) ./ abs.(us[:])) < 1e-9

    # Test 2D regular translation matrix
    V = regular_translation_matrix(medium, larger_order, order, ω, d)

    v1s = regular_basis_function(medium, ω, field_type)(larger_order,r)
    v2s = regular_basis_function(medium, ω, field_type)(order,r + d)

    @test maximum(abs.(V * v1s[:] - v2s[:]) ./ abs.(v2s[:])) < 1e-10

end

@testset "Translation 3D Matrices test" begin

    ω = rand() + 0.2
    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0)
    r = rand(3) - [0.5,0.5,0.5];
    d = rand(3) - [0.5,0.5,0.5];
    d = 10 * d * norm(r) / norm(d)

    # ω = 1.0;
    # r = [-0.16760003217875297,-0.35146260115553596,-0.2124509700641245];
    # d = [-2.9010760893418666, -3.353041498826335,0.12643092521032334];

    # Note that to be accurate the order of vs
    order = 4
    larger_order = 4*order
    basis_length = 2*larger_order+1
    
    # note the north pole is not well defined for spherical coordinates.
    basis = outgoing_basis_function(medium, ω, TractionType())
    @test_throws ArgumentError basis(order, [0.0,0.0,1.0])
    
    basis = regular_basis_function(medium, ω, DisplacementType())
    @test_throws ArgumentError basis(order, [0.0,0.0,1.0])

    # Test 3D outgoing translation matrix
    field_type = DisplacementType()

    V = regular_translation_matrix(medium, larger_order, order, ω, d, field_type)
    vs_in = regular_basis_function(medium, ω, field_type)(larger_order,r)
    vs_out = regular_basis_function(medium, ω, field_type)(order,r + d)

    @test maximum(abs.(V * transpose(vs_in) - transpose(vs_out))) < 1e-14
    @test norm(V * transpose(vs_in) - transpose(vs_out)) / norm(transpose(vs_out)) < 1e-14

    U = outgoing_translation_matrix(medium, larger_order, order, ω, d, field_type)
    vs_in = regular_basis_function(medium, ω, field_type)(larger_order,r)
    us_out = outgoing_basis_function(medium, ω, field_type)(order,r + d)

    @test maximum(abs.(U * transpose(vs_in) - transpose(us_out))) < 1e-13
    @test norm(U * transpose(vs_in) - transpose(us_out)) / norm(transpose(us_out)) < 1e-12
    
end