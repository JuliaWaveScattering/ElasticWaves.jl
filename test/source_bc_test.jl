using Test

@testset "Source matrices build test" begin

    ω = 1.0e6
    
    medium = Elastic(2; ρ = 7000.0, cp = 5000.0 - 0.0im, cs = 3500.0 - 0.0im)

    r1 = 1.0
    r2 = 1.5
    bearing = RollerBearing(medium = medium, 
        inner_radius = r1, outer_radius = r2
    )

    bc1 = TractionBoundary(inner=true)
    bc2 = TractionBoundary(outer=true)
    bc3 = DisplacementBoundary(outer = true)
    bc4 = DisplacementBoundary(inner = true)
    
    mode = rand(1:10)

    source_matrix1 = source_boundarycondition_system(ω, bearing, bc1, bc2, mode)
    source_matrix2 = source_boundarycondition_system(ω, bearing, bc3, bc2, mode)
    source_matrix3 = source_boundarycondition_system(ω, bearing, bc4, bc1, mode)
    source_matrix4 = source_boundarycondition_system(ω, bearing, bc4, bc3, mode)

    Minv_outer = boundarycondition_system(ω, bearing, bc3, bc2, mode)
    Minv_inner = boundarycondition_system(ω, bearing, bc4, bc1, mode)

    coes = ones(4)
    coes_inner = [1.0, 0.0, 1.0, 0.0]
    coes_outer = [0.0, 1.0, 0.0, 1.0]


    test_data1 = vcat(boundarycondition_mode(ω, bc1, bearing, mode)*coes_inner, boundarycondition_mode(ω, bc2, bearing, mode)*coes_outer)
    test_data2 = Minv_outer*coes_outer
    test_data3 = Minv_inner*coes_inner
    test_data4 = vcat(boundarycondition_mode(ω, bc4, bearing, mode)*coes_inner, boundarycondition_mode(ω, bc3, bearing, mode)*coes_outer)

    @test isapprox(norm(source_matrix1*coes - test_data1), 0.0; atol=1.0e-14)
    @test isapprox(norm(source_matrix2*coes - test_data2), 0.0; atol=1.0e-14)
    @test isapprox(norm(source_matrix3*coes - test_data3), 0.0; atol=1.0e-14)
    @test isapprox(norm(source_matrix4*coes - test_data4), 0.0; atol=1.0e-14)

end