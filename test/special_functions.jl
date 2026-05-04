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

    field_type = DisplacementType()
    # d = [2.0,2.0,1.0]
    # U = outgoing_translation_matrix(medium, 1, 1, 1.0, d, field_type);

    # l1 = 1; l2 =1;
    # len1 = (l1+1)^2
    # len2 = (l2+1)^2

    # dl = 1; dm = -1; l = 1; m = 1;
    # U[len1:(len1+2),len2:(len2+2)]


    ω = rand() + 0.2
    ω = 1.0;
    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0)
    r = rand(3) - [0.5,0.5,0.5];
    d = rand(3) - [0.5,0.5,0.5];
    d = 10 * d * norm(r) / norm(d)

    # Note that to be accurate the order of vs
    order = 4
    larger_order = 3*order
    basis_length = 2*larger_order+1
    
    # order = 2
    # larger_order = 2order
    # basis_length = larger_order+1
    
    # note the north pole is not well defined for spherical coordinates.
    basis = outgoing_basis_function(medium, ω, TractionType())
    @test_throws ArgumentError basis(order, [0.0,0.0,1.0])
    
    basis = regular_basis_function(medium, ω, DisplacementType())
    @test_throws ArgumentError basis(order, [0.0,0.0,1.0])

    #NOTE: we found that the regular basis in Mathematica is a bit different to the one in Julia. This needs checking next.
    # see test-addition-vector.jl
    # Test 3D outgoing translation matrix
    field_type = DisplacementType()
    V = regular_translation_matrix(medium, larger_order, order, ω, d, field_type)
    vs_in = regular_basis_function(medium, ω, field_type)(larger_order,r)
    vs_out = outgoing_basis_function(medium, ω, field_type)(order,r + d)

    V[end-2:end,1:3]

    # fs = rand(3 * (larger_order+1)^2)
    # vs * transpose(U) * fs
    # us * fs

    V * transpose(vs_in) - transpose(vs_out)
    

    field_type = PotentialType()
    U = outgoing_translation_matrix(medium, larger_order, order, ω, d, field_type)

    vs = regular_basis_function(medium, ω, field_type)(larger_order,r)
    us = outgoing_basis_function(medium, ω, field_type)(order,r + d)

    @test maximum(abs.(U * vs[:] - us[:]) ./ abs.(us[:])) < 1e-9

    # Test 3D regular translation matrix
    V = regular_translation_matrix(medium, larger_order, order, ω, d, field_type)

    v1s = regular_basis_function(medium, ω, field_type)(larger_order,r)
    v2s = regular_basis_function(medium, ω, field_type)(order,r + d)

    @test maximum(abs.(V * v1s[:] - v2s[:]) ./ abs.(v2s[:])) < 1e-10
end