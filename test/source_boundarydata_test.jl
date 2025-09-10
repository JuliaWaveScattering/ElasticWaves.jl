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

@testset "Source boundary_data test" begin

    medium=Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    #medium = Elastic(2; ρ = 1.0, cp = 1.5, cs = 1.0)
    r1 = 1.0
    r2 = 1.5
    bearing = RollerBearing(
        medium = medium,
        inner_radius = r1, outer_radius = r2
    )

    # we want wavelengths that go from about 1/10 of the bearing thickness radius to about 10 times the bearing radius
    krs = 0.1:0.1:10.0
    kps = krs .* (r1 + 0im)
    ωs = real.(kps .* medium.cp)
    kss = ωs ./ medium.cs

    # Source location
    θo = rand(0.0:0.01:2pi)
    ro = rand(bearing.inner_radius:0.01:bearing.outer_radius)
    xo, yo = radial_to_cartesian_coordinates([ro,θo])
    source_location = [xo,yo]
    p_map = SourceMap(source_location |> transpose, [1.0])
    s_map = SourceMap(source_location |> transpose, [0.0])

    # Defining modes to be used
    basis_order = 5
    sequecence = reshape(hcat([-(1:basis_order |> collect), 1:basis_order |> collect]...) |> transpose, 2*basis_order)
    modes = vcat([[0], sequecence]...)

    # Defining basic matrices to be used
    Mfor = [boundarycondition_system(ω, bearing, TractionBoundary(inner=true), TractionBoundary(outer=true), mode) for mode in modes, ω in ωs]
    Minv = [boundarycondition_system(ω, bearing, DisplacementBoundary(outer=true), TractionBoundary(outer=true), mode) for mode in modes, ω in ωs]

    Rt = [1.0 0.0 1.0 0.0;
        1.0 0.0 1.0 0.0;
        0.0 1.0 0.0 1.0;
        0.0 1.0 0.0 1.0]

    Ru = [0.0 1.0 0.0 1.0;
        0.0 1.0 0.0 1.0;
        0.0 1.0 0.0 1.0;
        0.0 1.0 0.0 1.0]

    Ts = [(Rt .* Mfor[i]) for i in CartesianIndices(Mfor)]
    Us = [(Ru .* Minv[i]) for i in CartesianIndices(Minv)]

    ## Checking formulas for the translation matrices
    p_medium = Acoustic(2; ρ = medium.ρ, c = medium.cp)
    s_medium = Acoustic(2; ρ = medium.ρ, c = medium.cs)

    n = basis_order
    maximum(abs.(regular_translation_matrix(p_medium, n, 0, ωs[1], [-xo,-yo])[1,:] .- [besselj.(i, kps[1]*ro) .* exp.(-1im*i*θo) for i in -n:n]))

    maximum(abs.(outgoing_translation_matrix(p_medium, n, 0, ωs[1], [-xo,-yo])[1,:] .- [hankelh1.(i, kps[1]*ro) .* exp.(-1im*i*θo) for i in -n:n]))


    U = [hankelh1(n,kps[i]*ro)*exp(-1im*n*θo) for n in -n:n, i in eachindex(ωs)]
    V = [besselj(n,kps[i]*ro)*exp(-1im*n*θo) for n in -n:n, i in eachindex(ωs)]

    U2=[outgoing_translation_matrix(p_medium,n,0,ωs[i],[-xo,-yo]) for i in eachindex(ωs) ]
    U2=vcat(U2...) |>transpose


    source_coeffs = (-im/4).* [hcat(U[:,i], V[:,i], zeros(length(modes)),  zeros(length(modes))) for i in eachindex(ωs)]

    tractions = [hcat([Ts[i,j]*source_coeffs[j][i,:] for i in eachindex(modes)]...) |> transpose for j in eachindex(ωs)]

    data_inv= [hcat([Us[i,j]*source_coeffs[j][i,:] for i in eachindex(modes)]...) |> transpose for j in eachindex(ωs)]
       
    #source_coeffs[1]
    # Recalculating tractions with new boundary data functions
    bd1 = [boundary_data(ωs[i], bearing, TractionBoundary(inner=true), p_map, s_map, modes).coefficients for i in eachindex(ωs)]
    bd2 = [boundary_data(ωs[i], bearing, TractionBoundary(outer=true), p_map, s_map, modes).coefficients for i in eachindex(ωs)]
    tractions_bd = [hcat(bd1[i],bd2[i]) for i in eachindex(ωs)]

    @test maximum([maximum(abs,(tractions_bd[i] .- tractions[i])/tractions[i]) for i in eachindex(ωs)]) <1e-9


    bd1_inv = [boundary_data(ωs[i], bearing, DisplacementBoundary(outer=true), p_map, s_map, modes).coefficients for i in eachindex(ωs)]
    bd2_inv = [boundary_data(ωs[i], bearing, TractionBoundary(outer=true), p_map, s_map, modes).coefficients for i in eachindex(ωs)]
    data_inv_bd = [hcat(bd1_inv[i],bd2_inv[i]) for i in eachindex(ωs)]

    @test maximum([maximum(abs,(data_inv[i] .- data_inv_bd[i])/data_inv[i]) for i in eachindex(ωs)]) <1e-9

end

@testset "Source boundary_data_test_no_dimension" begin


    medium=Elastic(2; ρ = 7800.0, cp = 5000.0, cs = 3500.0)
    #medium = Elastic(2; ρ = 1.0, cp = 1.5, cs = 1.0)
    r1 = 1.0
    r2 = 1.5
    bearing = RollerBearing(
        medium = medium,
        inner_radius = r1, outer_radius = r2
    )

    # we want wavelengths that go from about 1/10 of the bearing thickness radius to about 10 times the bearing radius
    krs = 0.1:0.1:10.0
    kps = krs .* (r1 + 0im)
    ωs = real.(kps .* medium.cp)
    kss = ωs ./ medium.cs
    
    non_dimensional_bearings=[nondimensionalise(bearing,ω) for ω in ωs]
    
    # Source location
    θo = rand(0.0:0.01:2pi)
    ro = rand(bearing.inner_radius:0.01:bearing.outer_radius)
    xo, yo = radial_to_cartesian_coordinates([ro,θo])
    source_location = [xo,yo]
    p_maps = [ SourceMap(source_location |> transpose,(ωs[i]/medium.cp)^2*[1.0]) for i in eachindex(ωs)]
    s_maps = [SourceMap(source_location |> transpose, (ωs[i]/medium.cp)^2*[0.0]) for i in eachindex(ωs)]

    # Defining modes to be used
    basis_order =5
    sequecence = reshape(hcat([-(1:basis_order |> collect), 1:basis_order |> collect]...) |> transpose, 2*basis_order)
    modes = vcat([[0], sequecence]...)

    

    # Defining basic matrices to be used
    Mfor = [boundarycondition_system(ωs[i], non_dimensional_bearings[i], TractionBoundary(inner=true), TractionBoundary(outer=true), mode) for mode in modes, i in eachindex(ωs)]
    Minv = [boundarycondition_system(ωs[i], non_dimensional_bearings[i], DisplacementBoundary(outer=true), TractionBoundary(outer=true), mode) for mode in modes, i in eachindex(ωs)]

    Rt = [1.0 0.0 1.0 0.0;
        1.0 0.0 1.0 0.0;
        0.0 1.0 0.0 1.0;
        0.0 1.0 0.0 1.0]

    Ru = [0.0 1.0 0.0 1.0;
        0.0 1.0 0.0 1.0;
        0.0 1.0 0.0 1.0;
        0.0 1.0 0.0 1.0]

    Ts = [(Rt .* Mfor[i]) for i in CartesianIndices(Mfor)]
    Us = [(Ru .* Minv[i]) for i in CartesianIndices(Minv)]

    ## Checking formulas for the translation matrices
    p_medium = Acoustic(2; ρ = medium.ρ, c = medium.cp)
    s_medium = Acoustic(2; ρ = medium.ρ, c = medium.cs)
    n = basis_order
    j=length(ωs)
    j=1

    p_medium_scaled = Acoustic(2; ρ = non_dimensional_bearings[j].medium.ρ, c = non_dimensional_bearings[j].medium.cp)
    s_medium_scaled = Acoustic(2; ρ = non_dimensional_bearings[j].medium.ρ, c = non_dimensional_bearings[j].medium.cs)

    maximum(abs.(hcat([regular_translation_matrix(p_medium, n, 0, ωs[j], [-xo,-yo])[1,:] .- [besselj.(i, kps[j]*ro) .* exp.(-1im*i*θo) for i in -n:n ] for j in eachindex(ωs) ]...)))

    
    #@test maximum(abs.(regular_translation_matrix(p_medium, n, 0, ωs[j], [-xo,-yo])[1,:] .- [besselj.(i, kps[j]*ro) .* exp.(-1im*i*θo) for i in -n:n ]) ) < 1e-7

    #@test maximum(abs.(outgoing_translation_matrix(p_medium, n, 0, ωs[j], [-xo,-yo])[1,:] .- [hankelh1.(i, kps[j]*ro) .* exp.(-1im*i*θo) for i in -n:n])) <1e-7
    
    @test maximum(abs.(regular_translation_matrix(p_medium_scaled, n, 0, ωs[j], [-xo,-yo])[1,:] .- [besselj.(i, ro) .* exp.(-1im*i*θo) for i in -n:n])) < 1e-7

    @test maximum(abs.(outgoing_translation_matrix(p_medium_scaled, n, 0, ωs[j], [-xo,-yo])[1,:] .- [hankelh1.(i, ro) .* exp.(-1im*i*θo) for i in -n:n])) <1e-7


    U = [hankelh1(n,ro)*exp(-1im*n*θo) for n in -n:n, i in eachindex(ωs)]
    V = [besselj(n,ro)*exp(-1im*n*θo) for n in -n:n, i in eachindex(ωs)]

    U2=[outgoing_translation_matrix(p_medium_scaled,n,0,ωs[i],[-xo,-yo]) for i in eachindex(ωs) ]
    U2=vcat(U2...) |>transpose


    source_coeffs = (-im/4).* [(ωs[i]/medium.cp)^2 .*hcat(U[:,i], V[:,i], zeros(length(modes)),  zeros(length(modes))) for i in eachindex(ωs)]

    tractions = [hcat([Ts[i,j]*source_coeffs[j][i,:] for i in eachindex(modes)]...) |> transpose for j in eachindex(ωs)]

    data_inv= [hcat([Us[i,j]*source_coeffs[j][i,:] for i in eachindex(modes)]...) |> transpose for j in eachindex(ωs)]
    
    #μ=medium.cs^2*medium.ρ
    #λ=medium.cp^2*medium.ρ-2*medium.cs^2*medium.ρ
    #source_coeffs[1]
    # Recalculating tractions with new boundary data functions
    bd1 = [boundary_data(ωs[i], non_dimensional_bearings[i], TractionBoundary(inner=true), p_maps[i], s_maps[i], modes).coefficients for i in eachindex(ωs)]
    bd2 = [boundary_data(ωs[i], non_dimensional_bearings[i], TractionBoundary(outer=true), p_maps[i], s_maps[i], modes).coefficients for i in eachindex(ωs)]
    tractions_bd = [hcat(bd1[i],bd2[i]) for i in eachindex(ωs)]
    #tractions_bd = [hcat(bd2[i],bd1[i]) for i in eachindex(ωs)]

    @test maximum([maximum(abs,tractions_bd[i] .- tractions[i]) for i in eachindex(ωs)]) <1e-7

    maximum(abs,tractions_bd[1])

    maximum(abs,tractions[1])

    bd1_inv = [ boundary_data(ωs[i], non_dimensional_bearings[i], DisplacementBoundary(outer=true), p_maps[i], s_maps[i], modes).coefficients for i in eachindex(ωs)]
    bd2_inv = [ boundary_data(ωs[i], non_dimensional_bearings[i], TractionBoundary(outer=true), p_maps[i], s_maps[i], modes).coefficients for i in eachindex(ωs)]
    data_inv_bd = [hcat(bd1_inv[i],bd2_inv[i]) for i in eachindex(ωs)]

    @test maximum([maximum(abs,data_inv[i] .- data_inv_bd[i]) for i in eachindex(ωs)]) <1e-7
end