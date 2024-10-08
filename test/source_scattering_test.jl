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

    field_reg = basis(basis_order, x - centre) * coes[:] 

    @test field_reg - field(source,x,ω) |> norm < 2e-14
end

@testset "Scattering" begin

    ω = 1.2
    basis_order = 9
    field_type = DisplacementType()
    order = basis_order
    T = Float64

    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)
    ks = ω / medium.cs
    kp = ω / medium.cp

    centre = [3.0, -3.0, 5.0]
    # centre = [0.0, 0.0, 0.0]

    source = plane_z_shear_source(medium)

    regular_coefficients = regular_spherical_coefficients(source)
    source_coes = regular_coefficients(basis_order,centre,ω)

    particle_medium = Elastic(3; ρ = 0.6, cp = 2.6, cs = 2.7 ./ 1.2)
    particle_shape = Sphere(centre,1.0)
    particle = Particle(particle_medium, particle_shape)

## Check scattered wave matches result from MultipleScattering

    Tmatrix = t_matrix(particle, medium, ω, basis_order)

    regular_coefficients = regular_spherical_coefficients(source)
    source_coes = regular_coefficients(basis_order,centre,ω)

    # coes_flat = (source_coes |> transpose)[:]
    coes_flat = (source_coes)[:]
    scattering_coes = Tmatrix * coes_flat
    # scattering_coes_reshape = (reshape(scattering_coes,3,:) |> transpose)[:]
    scattering_coes_reshape = scattering_coes

    outgoing_basis = outgoing_basis_function(medium, ω, DisplacementType())
    
    r = outer_radius(particle)
    x_vec =  [
        centre + spherical_to_cartesian_coordinates([2*r*rand() + r,pi * rand(),2pi * rand()])
    for i = 1:100 ]
    
    outside_fields = [outgoing_basis(basis_order, x - centre) * scattering_coes_reshape[:] for x in x_vec]
    
    sim = FrequencySimulation([particle],source)
    
    result = run(sim,x_vec,ω; basis_order = basis_order, only_scattered_waves = true)
    fs = field(result)
    
    a_vec = scattering_coes
    field_vec = field(sim, ω, x_vec, a_vec) - [field(sim.source)(x,ω) for x in x_vec]

    sim_a_vec = basis_coefficients(sim, ω; basis_order=basis_order)

    @test outside_fields - field_vec .|> norm |> maximum < 1e-14
    @test outside_fields - fs .|> norm |> maximum ≈ 0.0

## Check displacement boundary condition

    MGφΦs, MGχs = modal_system(particle, medium, ω, 1)
    MφΦ, GφΦ = MGφΦs[2]
    Mχ, Gχ = MGχs[2]

    function Tmat(l)
        TφΦ = (MGφΦs[l+1][1] \ MGφΦs[l+1][2])[[1,3],:]
        Tχ  = (MGχs[l+1][1] \ MGχs[l+1][2])[1]
        return [TφΦ zeros(T,2); T(0) T(0) Tχ]
    end

    function inner_mat(l)
        TφΦ = (MGφΦs[l+1][1] \ MGφΦs[l+1][2])[[2,4],:]
        Tχ  = (MGχs[l+1][1] \ MGχs[l+1][2])[2]
        return [TφΦ zeros(T,2); T(0) T(0) Tχ]
    end

    Tmatrix = t_matrix(particle, medium, ω, 1)
    @test Tmat(1) - Tmatrix[4:6,4:6] |> norm ≈ 0
    
    # The order l = 0 should be different
    @test Tmat(0) - Tmatrix[1:3,1:3] |> norm > 0.1

    in_matrix = internal_matrix(particle, medium, ω, 1)
    @test inner_mat(1) - in_matrix[4:6,4:6] |> norm ≈ 0

    # Test the internal_field
    gs = rand(3)
    bs = inner_mat(1) * gs 
    fs = Tmat(1) * gs 
    bs2 = inner_mat(1) * (Tmat(1) \ fs)
    @test bs - bs2 |> norm < 1e-12

    order = basis_order;

    pcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 
    Φcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 
    χcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 

    function sourceΦ_coes(order,centre,ω)
        return [pcoefs Φcoefs χcoefs] |> transpose
    end
    
    source_field = function (x1,ω) 
        source_basis = regular_basis_function(medium, ω, DisplacementType())
        source_basis(basis_order, x1 - centre) * sourceΦ_coes(basis_order,centre,ω)[:] 
    end
    
    sourceΦ = RegularSource{Elastic{3,T},WithoutSymmetry{3}}(medium, source_field, sourceΦ_coes)

    regular_coefficients = regular_spherical_coefficients(sourceΦ)
    source_coes = regular_coefficients(basis_order,centre,ω)

    # choose x on the boundary of the particle
    r = outer_radius(particle) + 10eps(T)
    xout = [
        centre + spherical_to_cartesian_coordinates([r, i * 2pi / 100, i * 7pi / 100]) 
    for i = 1:100] 
        
    r = outer_radius(particle) - 10eps(T)
    xin = [
        centre + spherical_to_cartesian_coordinates([r, i * 2pi / 100, i * 7pi / 100]) 
    for i = 1:100] 
    
    # xs = [centre + spherical_to_cartesian_coordinates([r,pi / 2,pi /2])] 

    # Need to include the source this time. 
    sim = FrequencySimulation([particle],sourceΦ)
    
    result = run(sim,xout,ω; basis_order = basis_order)
    fout = field(result)

    result = run(sim,xin,ω; basis_order = basis_order)
    fin = field(result)

    @test norm.(fin - fout) |> maximum < 1e-13

    # Alternative to calculate the field inside the particlce 
    in_matrix = internal_matrix(particle, medium, ω, basis_order)

    internal_coes = in_matrix * source_coes[:]

    basis = regular_basis_function(particle.medium, ω, DisplacementType())
    internal_fields = [basis(basis_order, x - centre) * internal_coes[:] for x in xout]

    @test norm.(internal_fields - fout) |> maximum < 1e-13

end

# check the traction boundary condition
@testset "Scattering traction" begin

    ω = 0.9
    basis_order = 10
    field_type = TractionType()

    order = basis_order
    T = Float64

    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)

    centre = [3.0, -3.0, 5.0]
    centre = [0.0, 0.0, 0.0]

    particle_medium = Elastic(3; ρ = 0.6, cp = 2.6, cs = 2.7 ./ 1.2)
    particle_shape = Sphere(centre,1.0)
    particle = Particle(particle_medium, particle_shape)

    order = basis_order;

    pcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 
    Φcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 
    χcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 

    # Φcoefs[1] = 0.0
    # χcoefs[1] = 0.0
    # pcoefs[1] = 0.0

    function sourceΦ_coes(order,centre,ω)
        return [pcoefs Φcoefs χcoefs] |> transpose
    end
    
    source_field = function (x1,ω) 
        source_basis = regular_basis_function(medium, ω, field_type)
        source_basis(basis_order, x1 - centre) * sourceΦ_coes(basis_order,centre,ω)[:] 
    end
    
    sourceΦ = RegularSource{Elastic{3,T},WithoutSymmetry{3}}(medium, source_field, sourceΦ_coes)

    regular_coefficients = regular_spherical_coefficients(sourceΦ)
    source_coes = regular_coefficients(basis_order,centre,ω)

    in_matrix = internal_matrix(particle, medium, ω, basis_order)
    internal_coes = in_matrix * source_coes[:]
    
    scat_matrix = t_matrix(particle, medium, ω, basis_order)
    external_coes = scat_matrix * source_coes[:]

    # choose x on the boundary of the particle
    r = outer_radius(particle)
    xs = [
        centre + spherical_to_cartesian_coordinates([r, i * 2pi / 100, i * 7pi / 100]) 
    for i = 1:100] 
         
    basis = regular_basis_function(particle.medium, ω, field_type)
    internal_fields = [basis(basis_order, x - centre) * internal_coes[:] for x in xs]
    
    basis = outgoing_basis_function(medium, ω, field_type)
    scat_fields = [basis(basis_order, x - centre) * external_coes[:] for x in xs]
    source_fields = [source_field(x,ω) for x in xs]
    
    external_fields = scat_fields + source_fields

    @test norm.(internal_fields - external_fields) |> maximum < 1e-13

end
