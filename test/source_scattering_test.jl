@testset "Source" begin

    ω = 1.2
    basis_order = 16
    field_type = DisplacementType()
    order = basis_order

    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)
    ks = ω / medium.cs

    source = plane_z_shear_source(medium)

    centre = [5.0, -3.0, 5.0]
    x = rand(3) .* 0.5 + centre
    
    regular_coefficients = regular_spherical_coefficients(source)
    coes = regular_coefficients(basis_order,centre,ω)
    
    basis = regular_basis_function(medium, ω, field_type)

    field_reg = basis(basis_order, x - centre) * coes[:] 

    @test norm(field_reg - field(source,x,ω)) / norm(field(source,x,ω)) < 1e-13

    # a monopoly source
    pos = [0.0, 0.0, 0.0]
    source =  regular_spherical_source(medium, rand(3); position = pos)

    regular_coefficients = regular_spherical_coefficients(source)
    coes = regular_coefficients(basis_order,centre,ω)
    
    basis = regular_basis_function(medium, ω, field_type)

    field_reg = basis(basis_order, x - centre) * coes[:] 

    @test norm(field_reg - field(source,x,ω)) / norm(field(source,x,ω)) < 1e-13
end

@testset "Single elastic particle scattering" begin

    ω = 1.2
    basis_order = 5
    order = basis_order
    field_type = DisplacementType()
    order = basis_order
    T = Float64

## comapre acoustic and elastic scattering for a pressure only source
    acoustic = Acoustic(3; ρ = 1.0, c = 1.0)
    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1e-8)
    ks = ω / medium.cs; kp = ω / medium.cp;

    centre = [0.0, 0.0, 0.0]
    particle_medium = Elastic(3; ρ = 0.1, cp = 0.4, cs = 1e-8)
    particle_shape = Sphere(centre,1.0)
    particle = Particle(particle_medium, particle_shape)
    
    particle_medium = Acoustic(3; ρ = 0.1, c = 0.4)
    particle_acoustic = Particle(particle_medium, particle_shape)

    # a pressure only source
    pos = [0.0, 0.0, 0.0]
    len = 3*(basis_order+1)^2
    # len = basisorder_to_basislength(Elastic{3}, basis_order)

    psource_coes = rand(Int(len/3));
    source_coes = [[c; 0.0; 0.0] for c in psource_coes];
    source_coes = vcat(source_coes...)

    Tmatrix = t_matrix(particle, medium, ω, basis_order)
    scat_coes = Tmatrix * source_coes
    
    Tmatrix = t_matrix(particle_acoustic, acoustic, ω, basis_order)
    scat_acoustic_coes = Tmatrix * psource_coes

    scat_acoustic_coes = [[c; 0.0; 0.0] for c in scat_acoustic_coes];
    scat_acoustic_coes = vcat(scat_acoustic_coes...)

    # convergence to a pure acoustic response is slow as cs -> 0, and leads to numerical instability. Essentially some shear waves get excited.
    @test norm(scat_acoustic_coes - scat_coes) / norm(scat_acoustic_coes) < 0.4

    # however the scattered acoustic waves are the same.
    @test norm(scat_acoustic_coes[1:3:end] - scat_coes[1:3:end]) / norm(scat_acoustic_coes[1:3:end]) < 5e-7
    @test norm(scat_acoustic_coes[1:3:12] - scat_coes[1:3:12]) / norm(scat_acoustic_coes[1:3:12]) < 5e-9

## test boundaru condition for general elastic source and particle
    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)
    ks = ω / medium.cs; kp = ω / medium.cp;

    centre = [3.0, -3.0, 5.0]

     # first using the source just to calculate appropriate scattering coefficients.
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

    particle_medium = Elastic(3; ρ = 0.6, cp = 2.6, cs = 2.7 ./ 1.2)
    particle_shape = Sphere(centre,1.0)
    particle = Particle(particle_medium, particle_shape)

## Check manually calculated scattered wave matches result from MultipleScattering

    Tmatrix = t_matrix(particle, medium, ω, basis_order)

    # coes_flat = (source_coes |> transpose)[:]
    coes_flat = (source_coes)[:]
    scattering_coes = Tmatrix * coes_flat

    outgoing_basis = outgoing_basis_function(medium, ω, field_type)
    
    r = outer_radius(particle)
    x_vec =  [
        centre + spherical_to_cartesian_coordinates([2*r*rand() + r,pi * rand(),2pi * rand()])
    for i = 1:100 ]
    
    outside_fields = [
        outgoing_basis(basis_order, x - centre) * scattering_coes 
    for x in x_vec]
    
    sim = FrequencySimulation([particle],sourceΦ)
    
    result = run(sim,x_vec,ω; basis_order = basis_order, only_scattered_waves = true)
    fs = field(result)
    
    a_vec = scattering_coes
    field_vec = field(sim, ω, x_vec, a_vec) - [field(sim.source)(x,ω) for x in x_vec]

    # sim_a_vec = basis_coefficients(sim, ω; basis_order=basis_order)

    @test outside_fields - field_vec .|> norm |> maximum < 1e-14
    @test outside_fields - fs .|> norm |> maximum ≈ 0.0

## Alternative to calculate the field inside particle without multiple scattering

    # choose x on the boundary of the particle
    # Need to avoid the north pole for spherical coordinates, as the basis functions are not well defined there.
        r = outer_radius(particle) + 10eps(T)
        xout = [
            centre + spherical_to_cartesian_coordinates([r, i * pi / 100, i * 7pi / 100]) 
        for i = 1:80] 
            
        r = outer_radius(particle) - 10eps(T)
        xin = [
            centre + spherical_to_cartesian_coordinates([r, i * pi / 100, i * 7pi / 100]) 
        for i = 1:80]

    particles = [particle] 
    sim = FrequencySimulation(particles,sourceΦ)

    result = run(sim,xout,ω; basis_order = basis_order)
    fout = field(result)

    in_matrix = internal_matrix(particle, medium, ω, basis_order)
    internal_coes = in_matrix * source_coes[:]

    basis = regular_basis_function(particle.medium, ω, field_type)
    internal_fields = [basis(basis_order, x - centre) * internal_coes[:] for x in xout]

    @test norm.(internal_fields - fout) |> maximum < 1e-13
end

@testset "Multiple elastic particle scattering" begin
    
    ω = 0.4
    basis_order = 6
    field_type = DisplacementType()
    order = basis_order
    T = Float64

    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)
    ks = ω / medium.cs; kp = ω / medium.cp;

    centre = [13.0, -3.0, 15.0]
    particle_medium = Elastic(3; ρ = 0.6, cp = 2.6, cs = 2.7 ./ 1.2)
    particle_shape = Sphere(centre,0.4)
    particle = Particle(particle_medium, particle_shape)

    centre2 = [0.0, 0.0, 0.0]
    particle_medium2 = Elastic(3; ρ = 1.6, cp = 1.6, cs = 0.9)
    particle_shape2 = Sphere(centre2,1.0)
    particle2 = Particle(particle_medium2, particle_shape2)
## Check displacement boundary condition

    # displacement source
        source = plane_z_shear_source(medium)
        regular_coefficients = regular_spherical_coefficients(source)
        source_coes = regular_coefficients(basis_order,centre,ω)

    # choose x on the boundary of the particle
    # need to avoid the North Pole (θ = 0) as the fields are not well defined there due to the singularity in the spherical coordinates.
        r = outer_radius(particle) + 10eps(T)
        xout = [
            centre + spherical_to_cartesian_coordinates([r, i * pi / 100, i * 7pi / 100]) 
        for i = 10:90]
            
        r = outer_radius(particle) - 10eps(T)
        xin = [
            centre + spherical_to_cartesian_coordinates([r, i * pi / 100, i * 7pi / 100]) 
        for i = 10:90]

        @test [norm(x - centre - [0.0, 0.0, 1.0]) for x in xout] |> minimum > 0.1
        @test [norm(x - centre - [0.0, 0.0, 1.0]) for x in xin] |> minimum > 0.1

        # The error in boundary will be dominated by the error made in the addition translation of the scattered waves from particle 2 to particle 1.
        U = outgoing_translation_matrix(medium, order, order, ω, origin(particle) - origin(particle2))

        errors = map(xout) do x
            vs = regular_basis_function(medium, ω, DisplacementType())(order, x - origin(particle))
            us = outgoing_basis_function(medium, ω, DisplacementType())(order, x - origin(particle2))
            norm(U * transpose(vs) - transpose(us)) / norm(transpose(us))
        end

        tol = maximum(errors)
        @test tol < 1e-8

        basis = regular_basis_function(medium, ω, field_type)
        errors = map(xout) do x
            source_reg = basis(basis_order, x - centre) * source_coes[:] 
            norm(source_reg - field(source,x,ω)) / norm(field(source,x,ω))
        end

        @test maximum(errors) < 1e-8
    
    # Need to include the source this time. Will also use two particles to check multiple scattering is working correctly.
        particles = [particle, particle2] 
        sim = FrequencySimulation(particles,source)

        result = run(sim,xout,ω; basis_order = basis_order)
        fout = field(result)

        result = run(sim,xin,ω; basis_order = basis_order)
        fin = field(result)

        # this test is currently failing:
        @test norm.(fin - fout) |> maximum < tol
        
    # An alternative way to calculate the field just outside the boundary of the particle is to use the scattering coefficients to calculate the exciting field, and then add the scattered field.

        f_arr = basis_coefficients(sim, ω; basis_order)
        fs = f_arr[:,1]

    # to recover the exciting field 
    # need to seperate the l=0 casef_arr = basis_coefficients(sim, ω; basis_order)
        t_mat = t_matrix(particle, medium, ω, order)
        L = length(fs)
        g0 = (t_mat[1,1] \ fs[1])

        # and now remove the l=0 cases
        # inds = [1, (order+1)^2 + 1, 2*(order+1)^2 + 1]
        inds = [1, 2, 3]

        exciting_coes = zeros(Complex{T}, L)
        exciting_coes[setdiff(1:end, inds)] = (t_mat[setdiff(1:end, inds), setdiff(1:end, inds)] \ fs[setdiff(1:end, inds)])
        exciting_coes[1] = g0

        # get scattering coefficients from exciting coes?
        @test abs.(fs - t_mat * exciting_coes) |> maximum < 1e-14

    reg_basis = regular_basis_function(medium, ω, field_type)
    field_reg = [reg_basis(order, x - centre) * exciting_coes[:] for x in xout]
    
    outgoing_basis = outgoing_basis_function(medium, ω, field_type)
    field_out = [outgoing_basis(order, x - centre) * fs for x in xout]

    fout2 = field_reg + field_out

    # this is how fin is calculated in the internal_field function, so naturally the test below should pass
    @test norm.(fout2 - fin) |> maximum < 1e-13

    # this test shows that the field calculated using the scattering coefficients, and the field calculated using the multiple scattering simulation are the same.

    basis = outgoing_basis_function(sim.source.medium, ω)
    fouts = map(eachindex(sim.particles)) do i
        p = sim.particles[i]
        [basis(order, x-origin(p)) * f_arr[:,i] for x in xout]
    end

    @test norm.(fouts[1] - field_out) |> maximum < 1e-16
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

    Φcoefs[1] = 0.0
    χcoefs[1] = 0.0

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
        centre + spherical_to_cartesian_coordinates([r, i * pi / 100, i * 7pi / 100]) 
    for i = 10:90]
         
    basis = regular_basis_function(particle.medium, ω, field_type)
    internal_fields = [basis(basis_order, x - centre) * internal_coes[:] for x in xs]
    
    basis = outgoing_basis_function(medium, ω, field_type)
    scat_fields = [basis(basis_order, x - centre) * external_coes[:] for x in xs]
    source_fields = [source_field(x,ω) for x in xs]
    
    external_fields = scat_fields + source_fields

    @test norm.(internal_fields - external_fields) |> maximum < 1e-13
    @test norm.(internal_fields - external_fields) ./ norm.(external_fields) |> maximum < 1e-13
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
