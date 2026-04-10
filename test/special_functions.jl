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
    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 0.5)
    r = rand(3) - [0.5,0.5,0.5];
    d = rand(3) - [0.5,0.5,0.5];
    d = 10 * d * norm(r) / norm(d)

    # Note that to be accurate the order of vs
    order = 4
    larger_order = 3*order
    basis_length = 2*larger_order+1
    
    # Test 3D outgoing translation matrix
    U = outgoing_translation_matrix(medium, larger_order, order, ω, d)

    field_type = PotentialType()
    vs = regular_basis_function(medium, ω, field_type)(larger_order,r)
    us = outgoing_basis_function(medium, ω, field_type)(order,r + d)

    @test maximum(abs.(U * vs[:] - us[:]) ./ abs.(us[:])) < 1e-9

    # Test 3D regular translation matrix
    V = regular_translation_matrix(medium, larger_order, order, ω, d)

    v1s = regular_basis_function(medium, ω, field_type)(larger_order,r)
    v2s = regular_basis_function(medium, ω, field_type)(order,r + d)

    @test maximum(abs.(V * v1s[:] - v2s[:]) ./ abs.(v2s[:])) < 1e-10
end