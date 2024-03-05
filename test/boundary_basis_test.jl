@testset "Boundary basis full test" begin

    medium = Elastic(2; ρ = 1.0, cp = 2.0, cs = 1.5)
    bearing = RollerBearing(medium = medium, inner_radius=1.0, outer_radius = 1.5)
    dr = bearing.outer_radius - bearing.inner_radius

    numberofsimulations = 6;

    # the higher the frequency, the worse the result. 
    ωs = LinRange(0.2,40,numberofsimulations) |> collect

    # the lower frequencies have less basis_order so can recover less information. For this reason we increase the numberofbasis with the frequency
    max_numberofbasis = 5
    numberofbasis = Int.(round.(LinRange(1,max_numberofbasis,numberofsimulations)))
    
    kps = (ωs ./ medium.cp)
    kps .* dr

## Forward problem with forcing on inner boundary to create a basis
    basis_order = 20;    
    θs = LinRange(0.0, 2pi, 2basis_order+2)[1:end-1]
    # using 2basis_order + 2 guarantees that we can exactly represent the above with Fourier modes with basis_order number of modes

    # choose a basis for the pressure and shear on the inner boundary
    # by using random data this is a hard test, because the fourier representation of random uniform data converges very slowly (because it is not smooth). This implies that all our methods converge slowly. Still, we should get precise results when disregarding fourier modes of higher orders that are not solved.

    random_basis() = [
        rand(length(θs)) .- 0.5 + (rand(length(θs)) .- 0.5) .* im
    for i = 1:max_numberofbasis]; 

    fp1s = random_basis(); fs1s = random_basis();
    fp2s = random_basis(); fs2s = random_basis();
    
## The forward problem

    # Create boundary data for the forward problem.
        bc1_forward = DisplacementBoundary(inner=true)
        bc2_forward = TractionBoundary(inner=true)

    # choose a combination of the boundary basis to be the boundary data

        boundarydata1s = map(eachindex(ωs)) do i
            fp1 = sum(fp1s[1:numberofbasis[i]]); 
            fs1 = sum(fs1s[1:numberofbasis[i]]);

            BoundaryData(bc1_forward, θs=θs, fields = hcat(fp1,fs1))
        end    

        boundarydata2s = map(eachindex(ωs)) do i
            fp2 = sum(fp2s[1:numberofbasis[i]]); 
            fs2 = sum(fs2s[1:numberofbasis[i]]);

            bd2_for = BoundaryData(bc2_forward, θs=θs, fields = hcat(fp2,fs2))
        end    

    # Solve the whole field for the forward problem        
        # the method specifies to use only stable modes.
        modal_method = ModalMethod(tol = 1e-11, only_stable_modes = true)

        forward_sims = map(eachindex(ωs)) do i
            BearingSimulation(ωs[i], bearing, boundarydata1s[i], boundarydata2s[i]; 
                method = modal_method,
                nondimensionalise = true
            )
        end;
        forward_waves = ElasticWave.(forward_sims);
        # [forward_waves[i].method.basis_order for i in eachindex(ωs)]

    # the lowest 3 frequencies have substantial errors below, as the problem become ill-posed for low basis_order and therefore does not completely describe the data given.
    forward_boundary_error = map(eachindex(ωs)) do i
        bds = boundary_data(forward_sims[i], forward_waves[i]);
        error1 = norm(boundarydata1s[i].fields - bds[1].fields) / norm(boundarydata1s[i].fields);
        error2 = norm(boundarydata2s[i].fields - bds[2].fields) / norm(boundarydata2s[i].fields);

        max(error1,error2)
    end

    # if we eliminate the basis_orders that the forward problem could not solve for, then the recovery is well below the specified tolerance (method.tol) as it should be
    forward_boundary_error = map(eachindex(ωs)) do i
        bd1 = fields_to_fouriermodes(boundarydata1s[i], forward_waves[i].method.basis_order)
        bd1 = fouriermodes_to_fields(bd1);
        
        bd2 = fields_to_fouriermodes(boundarydata2s[i], forward_waves[i].method.basis_order)
        bd2 = fouriermodes_to_fields(bd2);

        bds = boundary_data(forward_sims[i], forward_waves[i]);

        error1 = norm(bd1.fields - bds[1].fields) / norm(bd1.fields)
        error2 = norm(bd2.fields - bds[2].fields) / norm(bd2.fields)
        max(error1,error2)
    end
    @test maximum(forward_boundary_error) < modal_method.tol

## The inverse problem
        bc1_inv = DisplacementBoundary(outer=true)
        bc2_inv = TractionBoundary(outer=true)

    # Each sensor makes two measurements. For this test we want exact recovery. So we choose to have at least as many measurements as unknowns. The number of unknowns is equal to the number of basis.
        numberofsensors = numberofbasis
        
        θs_inv = [
            LinRange(0, 2pi, numberofsensors[i] + 1)[1:end-1] 
        for i in eachindex(ωs)]

        boundarydata1_inverses = [ 
            BoundaryData(bc1_inv, bearing.outer_radius, θs_inv[i], forward_waves[i])
        for i in eachindex(ωs)];

        boundarydata2_inverses = [ 
            BoundaryData(bc2_inv, bearing.outer_radius, θs_inv[i], forward_waves[i])
        for i in eachindex(ωs)];

    # will use a fine grid to demonstrate and explain some issues.    
        boundarydata1_fine_inverses = [ 
            BoundaryData(bc1_inv, bearing.outer_radius, θs, forward_waves[i])
        for i in eachindex(ωs)];

        boundarydata2_fine_inverses = [ 
            BoundaryData(bc2_inv, bearing.outer_radius, θs, forward_waves[i])
        for i in eachindex(ωs)];

    # Create a boundary basis for the inverse problem for the inner boundaries
        boundarybasis1_inverses = map(eachindex(ωs)) do i
            bd1s = [
                BoundaryData(bc1_forward, θs = θs, fields = hcat(fp1s[j],fs1s[j]))
            for j in 1:numberofsensors[i]];
            BoundaryBasis(bd1s)
        end;

        boundarybasis2_inverses = map(eachindex(ωs)) do i
            bd1s = [
                BoundaryData(bc2_forward, θs = θs, fields = hcat(fp1s[j],fs1s[j]))
            for j in 1:numberofsensors[i]];
            BoundaryBasis(bd1s)
        end;

## Inverse prior methods 
    # we solve the inverse problem in three different ways: 
    # - method 1 : with no prior knowledge,
    # - method 2 : with prior for one boundary condition
    # - method 3 : with prior for two boundary conditions

    # To test the methods we need the true data from the forward modals. Better to use data from the forward waves, rather than the data given to the forward problem, as discussed before the forward problem was not able to resolve all the modes of the data.
    forward_inner_data = [
        boundary_data(forward_sims[i], forward_waves[i])
    for i in eachindex(ωs)];

## Method 1
    # when using very limited data, this method performs poorly. It requires fine data to get fine predictions
    
    # using coarse boundary data boundarydata1_inverses for the modal method will lead to bad results, as there is just not enough data to calculate even the modes that were accurately resolved by the forward problem
        # sims = map(eachindex(ωs)) do i
        #     BearingSimulation(ωs[i], modal_method, bearing, boundarydata1_inverses[i], boundarydata2_inverses[i];
        #         nondimensionalise = true
        #     )
        # end;
        # inverse_waves = ElasticWave.(sims);

    # so instead we run the modal method with a lot of data just to demonstrate that it does work
    sims = map(eachindex(ωs)) do i
        BearingSimulation(ωs[i], modal_method, bearing, boundarydata1_fine_inverses[i], boundarydata2_fine_inverses[i];
            nondimensionalise = true
        )
    end;
    inverse_waves = ElasticWave.(sims);
    
    predicted_inner_data = [
        boundary_data(forward_sims[i], inverse_waves[i])
    for i in eachindex(ωs)];

    inner_boundary_errors = map(eachindex(ωs)) do i
        error1 = norm(forward_inner_data[i][1].fields - predicted_inner_data[i][1].fields) / norm(forward_inner_data[i][1].fields)
        error2 = norm(forward_inner_data[i][2].fields - predicted_inner_data[i][2].fields) / norm(forward_inner_data[i][2].fields)
        max(error1,error2)
    end

    @test maximum(inner_boundary_errors) < 2e-8


## Method 2

    # here we use coarse boundary data for the boundary with the prior. For the other boundary we need to use fine boundary data, otherwise the prior method is not even an improvement on the ModalMethod
    modal_method = ModalMethod(tol = 1e-11, only_stable_modes = true)
    method = PriorMethod(tol = modal_method.tol, modal_method = modal_method)

    sims = map(eachindex(ωs)) do i
        BearingSimulation(ωs[i], method, bearing, 
            boundarydata1_inverses[i], 
            boundarydata2_fine_inverses[i];
            boundarybasis1 = boundarybasis1_inverses[i]
        )
    end;
    inverse_waves = ElasticWave.(sims);

    # Check the error on the boundaries of the forward problem
    inner_boundary_errors = map(eachindex(ωs)) do i
        predict_bds = boundary_data(forward_sims[i], inverse_waves[i])

        error1 = norm(forward_inner_data[i][1].fields - predict_bds[1].fields) / norm(forward_inner_data[i][1].fields)
        error2 = norm(forward_inner_data[i][2].fields - predict_bds[2].fields) / norm(forward_inner_data[i][2].fields)
        max(error1,error2)
    end

    # method should always work for high enough frequencies
    @test inner_boundary_errors[end] < 1e-10

    # There error is surprisingly large, I find. What has happened is that Method 2 was solved to a higher basis_order (more resolution) than the forward problem. This means this Method 2 used parts of the boundarybasis1_inverses that the forward problem ignored (disgarded). So, for some frequencies, Method 2 is closer to the original forward boundary data than it is to the inner boundary data predicted by the forward problem
    original_inner_boundary_errors = map(eachindex(ωs)) do i
        predict_bds = boundary_data(forward_sims[i], inverse_waves[i])

        error1 = norm(boundarydata1s[i].fields - predict_bds[1].fields) / norm(boundarydata1s[i].fields)
        error2 = norm(boundarydata2s[i].fields - predict_bds[2].fields) / norm(boundarydata2s[i].fields)

        max(error1,error2)
    end
    @test original_inner_boundary_errors[end] < 1e-10
    
    forward_inner_boundary_errors = map(eachindex(ωs)) do i

        error1 = norm(boundarydata1s[i].fields - forward_inner_data[i][1].fields) / norm(boundarydata1s[i].fields)
        error2 = norm(boundarydata2s[i].fields - forward_inner_data[i][2].fields) / norm(boundarydata2s[i].fields)

        max(error1,error2)
    end

    # There is no way to perfectly recover the original inner boundary data for all frequencies, as the forward problem is not perfectly solved for all frequencies. However, we can exactly recover the solution of the forward problem by solving Method 2 with the same basis_order as the forward problem

    sims = map(eachindex(ωs)) do i
        modal_method = ModalMethod(tol = 1e-11, only_stable_modes = true, basis_order = forward_waves[i].method.basis_order)
        method = PriorMethod(tol = modal_method.tol, modal_method = modal_method)

        BearingSimulation(ωs[i], method, bearing, 
            boundarydata1_inverses[i], 
            boundarydata2_fine_inverses[i];
            boundarybasis1 = boundarybasis1_inverses[i]
        )
    end;
    inverse_waves = ElasticWave.(sims);

    inner_boundary_errors = map(eachindex(ωs)) do i
        predict_bds = boundary_data(forward_sims[i], inverse_waves[i])

        error1 = norm(forward_inner_data[i][1].fields - predict_bds[1].fields) / norm(forward_inner_data[i][1].fields)
        error2 = norm(forward_inner_data[i][2].fields - predict_bds[2].fields) / norm(forward_inner_data[i][2].fields)
        max(error1,error2)
    end

    @test maximum(inner_boundary_errors) < 1e-10


## Method 3

    sims = map(eachindex(ωs)) do i
        # modal_method = ModalMethod(tol = 1e-11, only_stable_modes = true)
        modal_method = ModalMethod(tol = 1e-11, only_stable_modes = true, basis_order = forward_waves[i].method.basis_order)
        method = PriorMethod(tol = modal_method.tol, modal_method = modal_method)

        BearingSimulation(ωs[i], method, bearing, 
            boundarydata1_inverses[i], 
            boundarydata2_inverses[i];
            boundarybasis1 = boundarybasis1_inverses[i],
            boundarybasis2 = boundarybasis2_inverses[i],
        )
    end;
    inverse_waves = ElasticWave.(sims);

    map(eachindex(ωs)) do i
        (ωs[i], forward_waves[i].method.basis_order)
    end

    # inverse_sim = BearingSimulation(ω, method, bearing, bd1_inverse, bd2_inverse;
    #     boundarybasis1 = boundarybasis1,
    #     nondimensionalise = true
    # );

    # inverse_wave = ElasticWave(inverse_sim);

# ## Test how well we recover the inner traction
#     x_inner = [
#         radial_to_cartesian_coordinates([bearing.inner_radius,θ]) 
#     for θ in θs]
    
#     field1_inner = [
#         field(inverse_wave, x, bc1_forward.fieldtype) 
#     for x in x_inner]
#     field1_inner = hcat(field1_inner...) |> transpose |> collect


#     # even excluding some modes the recovery is very good with just one sensor
#     # this is because the higher order modes (in this case) contribute very little to the field
#     @test norm(field1_inner - bd1_for.fields) / norm(bd1_for.fields) < 1e-8
    
#     @test norm(wave.potentials[1].coefficients - inverse_wave.potentials[1].coefficients) / norm(wave.potentials[1].coefficients) < 2e-12

#     @test norm(wave.potentials[2].coefficients - inverse_wave.potentials[2].coefficients) / norm(wave.potentials[2].coefficients) < 2e-12

#     modes1 = field_modes(wave, bearing.inner_radius, bc1_forward.fieldtype);
#     inverse_modes1 = field_modes(inverse_wave, bearing.inner_radius, bc1_forward.fieldtype);

#     @test norm(modes1 - inverse_modes1) / norm(modes1) < 2e-12
end


# inner_boundary_errors = map(eachindex(ωs)) do i
#     # order = inverse_waves[i].method.modal_method.basis_order
#     order = forward_waves[i].method.basis_order

#     bd1 = fields_to_fouriermodes(forward_inner_data[i][1], order)
#     bd1 = fouriermodes_to_fields(bd1);
    
#     bd2 = fields_to_fouriermodes(forward_inner_data[i][2], order)
#     bd2 = fouriermodes_to_fields(bd2);

#     predict_bds = boundary_data(forward_sims[i], inverse_waves[i])
    
#     predict_bd1 = fields_to_fouriermodes(predict_bds[1], order)
#     predict_bd1 = fouriermodes_to_fields(predict_bd1);

#     predict_bd2 = fields_to_fouriermodes(predict_bds[2], order)
#     predict_bd2 = fouriermodes_to_fields(predict_bd2);

#     error1 = norm(bd1.fields - predict_bds[1].fields) / norm(bd1.fields)
#     error2 = norm(bd2.fields - predict_bds[2].fields) / norm(bd2.fields)

#     max(error1,error2)
# end