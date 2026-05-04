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

    @test field_reg - field(source,x,ω) |> norm < 2e-14
end

@testset "Single elastic particle scattering" begin

    ω = 1.2
    basis_order = 6
    order = basis_order
    field_type = DisplacementType()
    order = basis_order
    T = Float64

    medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 1.0 ./ 1.2)
    ks = ω / medium.cs; kp = ω / medium.cp;

    centre = [3.0, -3.0, 5.0]

     # first using the source just to calculate appropriate scattering coefficients.
    pcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 
    Φcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 
    χcoefs = [rand() + rand()*im for l = 0:order for m = -l:l] 

    function sourceΦ_coes(order,centre,ω)
        return [pcoefs; Φcoefs; χcoefs]
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

## given the scattering coefficients, can we recover the exciting field and internal field correctly?

    # first create appropriate scattering coefficients
    t_mat = t_matrix(particle, medium, ω, order)
    source_coes = regular_coefficients(order,centre,ω)
    scat_coes = t_mat * source_coes[:]
    
    # to recover the exciting field 
    # need to seperate the l=0 case
        # fs = deepcopy(scat_coes)
        # L = length(fs)
        # g0 = (t_mat[1,1] \ fs[1])

        # # and now remove the l=0 cases
        # inds = [1, (order+1)^2 + 1, 2*(order+1)^2 + 1]
        # t_mat = t_mat[setdiff(1:end, inds), setdiff(1:end, inds)]

        # exciting_coefs = zeros(Complex{T}, L)
        # exciting_coefs[setdiff(1:end, inds)] = (t_mat \ fs[setdiff(1:end, inds)])
        # exciting_coefs[1] = g0

        # # check are the same as the source coefficients
        # # note we do not check the l=0 coefficients for the shear waves, as these have to be zero.
        # @test abs.(exciting_coefs[setdiff(1:end, inds[2:end])] - source_coes[setdiff(1:end, inds[2:end])]) |> maximum < 1e-10
    
    # given scat_coes, can we recover the exciting field and internal field correctly?
    # fin = [internal_field(x, particle, sourceΦ, ω, scat_coes, field_type) for x in xin]
    
    # # put it all together and check boundary condition: fout == fin?
    # reg_basis = regular_basis_function(medium, ω, field_type)
    # field_reg = [reg_basis(order, x - centre) * exciting_coefs[:] for x in xout]
    
    # outgoing_basis = outgoing_basis_function(medium, ω, field_type)
    # field_out = [outgoing_basis(order, x - centre) * scat_coes[:] for x in xout]

    # fout = field_reg + field_out

    # @test norm.(fout - fin) |> maximum < 1e-13
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

    # Test the internal_field
    # gs = rand(3)
    # bs = inner_mat(1) * gs 
    # fs = Tmat(1) * gs 
    # bs2 = inner_mat(1) * (Tmat(1) \ fs)
    # @test bs - bs2 |> norm < 1e-12

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
            vs = regular_basis_function(medium, ω, PotentialType())(order, x - origin(particle))
            us = outgoing_basis_function(medium, ω, PotentialType())(order, x - origin(particle2))
            norm(U * vs[:] - us[:]) / norm(us[:])
        end

        tol = maximum(errors)
        @test tol < 1e-9

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
        norm.(fin)
        norm.(fout)

        fin[end] .|> abs
        fout[end] .|> abs


    # An alternative way to calculate the field just outside the boundary of the particle is to use the scattering coefficients to calculate the exciting field, and then add the scattered field.

        f_arr = basis_coefficients(sim, ω; basis_order)
        fs = f_arr[:,1]

    # to recover the exciting field 
    # need to seperate the l=0 casef_arr = basis_coefficients(sim, ω; basis_order)
        t_mat = t_matrix(particle, medium, ω, order)
        L = length(fs)
        g0 = (t_mat[1,1] \ fs[1])

        # and now remove the l=0 cases
        inds = [1, (order+1)^2 + 1, 2*(order+1)^2 + 1]

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

    # fouts[2] + source == field_reg ?
    # Are these the same coefficients?
    
    U = outgoing_translation_matrix(medium, order, order, ω, origin(particles[1]) - origin(particles[2]))
    reshape(f_arr[:,2],:,3)
    exciting_coes2 = transpose(U) * f_arr[:,2] + source_coes[:]
    exciting_coes2 = reshape(exciting_coes2, :, 3)
    exciting_coes = reshape(exciting_coes, :, 3)
    
    # doesnt match due to the l=0 shear coefficients
    abs.(exciting_coes - exciting_coes2)[1,1:3] 
    
    # however if we remove the l = 0 shear coefficients, which make no contribution to the field, then the remaining coefficients are the same.
    @test abs.(exciting_coes[2:end,[1,2,3]] - exciting_coes2[2:end,[1,2,3]]) |> maximum < 1e-10
    @test abs.(exciting_coes[1,1] - exciting_coes2[1,1]) < 1e-14

    @test abs.(t_mat * exciting_coes[:] - t_mat * exciting_coes2[:]) |> maximum < 1e-14

    field_reg = [reg_basis(order, x - centre) * exciting_coes[:] for x in xout]
    field_reg2 = [reg_basis(order, x - centre) * exciting_coes2[:] for x in xout]

    norm.(field_reg - field_reg2) |> maximum

    out_order = 1
    out_order = order

    field_type2 = PotentialType()
    field_type2 = DisplacementType()

    U = outgoing_translation_matrix(medium, out_order, order, ω, origin(particles[1]) - origin(particles[2]))
    exciting_no_source_coes = transpose(U) * f_arr[:,2]

    reg_basis = regular_basis_function(medium, ω, field_type2)
    exciting_no_source = [reg_basis(out_order, x - centre) * exciting_no_source_coes[:] for x in xout]
   
    basis = outgoing_basis_function(sim.source.medium, ω, field_type2)
    
    f2_out = reshape(f_arr[:,2],:,3)[1:(out_order+1)^2,:]
    exciting_no_source_exact = map(xout) do x
        basis(out_order, x-origin(particle2)) * f2_out[:]
    end

    x = xout[80]
    data1 = reg_basis(out_order, x - centre) * transpose(U) 
    data2 = basis(out_order, x-origin(particle2))

    data = abs.(data1 - data2)
    data = reshape(data, 3,(out_order+1)^2, 3)
    out_basis = reshape(data2, 3, (out_order+1)^2, 3)

    pdata = data[:,:,1]
    φdata = data[:,:,2]
    χdata = data[:,:,3]
    norm(data) / norm(out_basis)

    norm.(exciting_no_source - exciting_no_source_exact) ./ norm.(exciting_no_source_exact)  

    norm.(exciting_no_source_exact)
    norm.(exciting_no_source)  

    data1n = [data1[:,(out_order+1)^2 + 1:((out_order+1)^2*2)][:,i] |> norm for i in 1:((out_order+1)^2)]
    data2n = [data2[:,(out_order+1)^2 + 1:((out_order+1)^2*2)][:,i] |> norm for i in 1:((out_order+1)^2)]
    
    outgoing_basis = outgoing_basis_function(medium, ω, field_type)
    field_out = [outgoing_basis(order, x - centre) * fs for x in xout]

    fout2 = field_reg + field_out

    @test norm.(fout2 - fout) |> maximum < 1e-13

    in_matrix = internal_matrix(particle, medium, ω, basis_order)
    internal_coes = in_matrix * exciting_coes2[:]
    abs.(internal_coes - in_matrix * exciting_coes[:]) |> maximum

    reshape(internal_coes, :, 3)
    
    # note the non-zero l=0 coefficients for the shear wave are not the same.
    reshape(U * f_arr[:,2], :, 3)

    abs.(t_mat * exciting_coes[:] - t_mat * exciting_coes2[:]) |> maximum



    

    t_matrices = get_t_matrices(sim.source.medium, sim.particles, ω, basis_order)
    S = scattering_matrix(sim.source.medium, sim.particles, t_matrices, ω, basis_order)

    source_coefficient = regular_spherical_coefficients(sim.source)
    forcing = reduce(vcat, [source_coefficient(basis_order,origin(p),ω)[:] for p in sim.particles])

    # Find scattering coefficients by solving this forcing
    a = (S + I) \ forcing
    (S + I)*a - forcing |> norm

    a = reshape(a,:,length(sim.particles))
    f_arr2 = deepcopy(a)
    for i in axes(a,2)
        f_arr2[:,i] = t_matrices[i] * a[:,i]
    end
    
    f_arr = basis_coefficients(sim, ω; basis_order)
    f_arr2 - f_arr

    function S_block(j,l)
        if j == l
            return sparse(zeros(Complex{T}, N, N))
        else
            x_lj = origin(particles[j]) .- origin(particles[l])
            U = outgoing_translation_matrix(medium, basis_order, basis_order, ω, x_lj)
            return - transpose(U) * t_matrices[l]
        end
    end

    U = outgoing_translation_matrix(medium, order, order, ω, origin(particles[1]) - origin(particles[2]))
    exciting_coes2 = transpose(U) * f_arr[:,2]

    maximum(abs, reshape(S*a,:,2)[:,1] - S_block(1,2) * reshape(a,:,2)[:,2])

    maximum(abs, reshape(S*a,:,2)[:,1] + exciting_coes2)

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
        return [pcoefs; Φcoefs; χcoefs]
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
