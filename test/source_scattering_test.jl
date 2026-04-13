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
    ω = 0.01
    basis_order = 9
    basis_order = 6
    field_type = DisplacementType()
    order = basis_order
    T = Float64

    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)
    ks = ω / medium.cs
    kp = ω / medium.cp

    centre = [3.0, -3.0, 5.0]
    source = plane_z_shear_source(medium)

    regular_coefficients = regular_spherical_coefficients(source)
    source_coes = regular_coefficients(basis_order,centre,ω)

    particle_medium = Elastic(3; ρ = 0.6, cp = 2.6, cs = 2.7 ./ 1.2)
    particle_shape = Sphere(centre,1.0)
    particle = Particle(particle_medium, particle_shape)

    particle_medium2 = Elastic(3; ρ = 1.6, cp = 1.6, cs = 0.9)
    centre2 = [0.0, 0.0, 0.0]
    particle_shape2 = Sphere(centre2,0.7)
    particle2 = Particle(particle_medium2, particle_shape2)

## Check scattered wave matches result from MultipleScattering

    Tmatrix = t_matrix(particle, medium, ω, basis_order)

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
    
    outside_fields = [
        outgoing_basis(basis_order, x - centre) * scattering_coes_reshape[:] 
    for x in x_vec]
    
    sim = FrequencySimulation([particle],source)
    
    result = run(sim,x_vec,ω; basis_order = basis_order, only_scattered_waves = true)
    fs = field(result)
    
    a_vec = scattering_coes
    field_vec = field(sim, ω, x_vec, a_vec) - [field(sim.source)(x,ω) for x in x_vec]

    sim_a_vec = basis_coefficients(sim, ω; basis_order=basis_order)

    @test outside_fields - field_vec .|> norm |> maximum < 1e-14
    @test outside_fields - fs .|> norm |> maximum ≈ 0.0

## Check displacement boundary condition

    # Test the internal_field
    # gs = rand(3)
    # bs = inner_mat(1) * gs 
    # fs = Tmat(1) * gs 
    # bs2 = inner_mat(1) * (Tmat(1) \ fs)
    # @test bs - bs2 |> norm < 1e-12

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
    
    # Need to include the source this time. Will also use two particles to check multiple scattering is working correctly.
    particles = [particle, particle2] 
    sim = FrequencySimulation(particles,sourceΦ)

    # U = outgoing_translation_matrix(medium, basis_order, basis_order, ω, xout[1])
    # transpose(U) * Tmatrix

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

    centre = [0.0, 0.0, 0.0]
    particle_medium = Elastic(3; ρ = 0.6, cp = 2.6, cs = 2.7 ./ 1.2)
    particle_shape = Sphere(centre,1.0)
    particle = Particle(particle_medium, particle_shape)
    
    particle_medium2 = Elastic(3; ρ = 1.6, cp = 1.6, cs = 0.9)
    centre = [3.0, -3.0, 5.0]
    particle_shape2 = Sphere(centre,0.7)
    particle2 = Particle(particle_medium2, particle_shape2)

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

@testset "2D T-matrix" begin
    ω = 1.2
    basis_order = 9
    field_type = TractionType()
    order = basis_order
    T = Float64

    medium = Elastic(2; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)
    ks = ω / medium.cs
    kp = ω / medium.cp

    #centre = [3.0, -3.0, 5.0]
    centre = [0.0, 0.0]

    particle_medium = Elastic(2; ρ = 0.6, cp = 2.6, cs = 2.7 ./ 1.2)
    particle_shape = Sphere(centre,1.0)
    particle = Particle(particle_medium, particle_shape)

    T = t_matrix(particle, medium, ω, basis_order)
    modes = -basis_order:basis_order

    # The incident field coefficients for a plane wave in the x direction is given by i^n, where n is the mode number. 
    incident_coefficients = [[ComplexF64(im)^n, ComplexF64(im)^n] for n in modes]
    incident_coefficients = vcat(incident_coefficients...)

    scattered_coefficients = T * incident_coefficients

    m13, m24 = modal_system(particle, medium, ω, basis_order)

    @test maximum(abs.(m13 * incident_coefficients + m24 * scattered_coefficients)) < 1e-13
end
