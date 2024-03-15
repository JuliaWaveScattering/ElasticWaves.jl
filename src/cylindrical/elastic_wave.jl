function boundarycondition_mode(ω::AbstractFloat, bc::BoundaryCondition, bearing::RollerBearing, mode::Int)
    r = (bc.inner == true) ? bearing.inner_radius : bearing.outer_radius
    return hcat(
        pressure_field_mode(ω, r, bearing.medium, mode, bc.fieldtype),
        shear_field_mode(ω, r, bearing.medium, mode, bc.fieldtype)
    )
end

function boundarycondition_system(ω::AbstractFloat, bearing::RollerBearing, bc1::BoundaryCondition, bc2::BoundaryCondition, mode::Int)

    if bc1 == bc2
        error("Both boundary conditions are the same. It is not possible to solve the system with just one type of boundary data/condition.")
    end     

    first_boundary = boundarycondition_mode(ω, bc1, bearing, mode)
    second_boundary = boundarycondition_mode(ω, bc2, bearing, mode)

    return vcat(first_boundary, second_boundary)
end

function ElasticWave(sim::BearingSimulation)

    # using the bearing properties before nondimensionalise
    bearing = sim.bearing;
    simcopy = deepcopy(sim);

    if simcopy.nondimensionalise
        nondimensionalise!(simcopy)
    end;

    ω = simcopy.ω;

    modes, coefficients = modes_coefficients!(simcopy);

    method = simcopy.method

    kP = ω / bearing.medium.cp;
    kS = ω / bearing.medium.cs;

    coefficients = if simcopy.nondimensionalise
        coefficients ./ kP^2
    else coefficients   
    end

    pressure_coefficients = coefficients[:,1:2] |> transpose
    shear_coefficients = coefficients[:,3:4] |> transpose

    φ = HelmholtzPotential{2}(bearing.medium.cp, kP, pressure_coefficients, modes)
    ψ = HelmholtzPotential{2}(bearing.medium.cs, kS, shear_coefficients, modes)

    return ElasticWave(ω, bearing.medium, [φ, ψ], method)
end

function modes_coefficients!(sim::BearingSimulation{ModalMethod})

    ω = sim.ω
    bearing = sim.bearing
    method = sim.method

    T = typeof(ω)

    modes = sim.method.modes

    coefficients = [zeros(Complex{T}, 4) for n = modes]

    mode_errors = zeros(T, length(modes))

    is = sortperm(modes, by = abs)

    for i in is
        A = boundarycondition_system(ω, bearing, sim.boundarydata1.boundarytype, sim.boundarydata2.boundarytype, modes[i])
        b = [
            sim.boundarydata1.coefficients[i,:];
            sim.boundarydata2.coefficients[i,:]
        ]

        # solve A*x = b with tikinov regulariser
        # x = A \ b
        δ = method.regularisation_parameter
        x = [A; sqrt(δ) * I] \ [b; zeros(size(A)[2])]

        relative_error_x = cond(A) * eps(T) 
        
        error = if norm(b) > 0
            relative_error = norm(A*x - b) / norm(b);
            max(relative_error,relative_error_x)
        else relative_error_x
        end

        coefficients[i] = x
        mode_errors[i] = error

        if error > method.tol
            @warn "The relative error for the boundary conditions was $(error) for (ω,mode) = $((ω,modes[i]))"
            if method.only_stable_modes 
                mode_errors[i] = zero(T)
                break 
            end
        end
    end

    coefficients = transpose(hcat(coefficients...)) |> collect

    # just in case the loop terminated early due to keeping only_stable_modes
    i = findfirst(mode_errors[is] .== zero(T))
    if !isnothing(i)
        if i == 1
            error("there are no stable modes, and therefore no solution found")
        end    
        
        # best to keep the indices in the same order they were given and not in the sorted order
        inds = sort(is[1:i-1])

        coefficients = coefficients[inds,:]
        mode_errors = mode_errors[inds]
        modes = modes[inds]
    end       
    
    @reset method.mode_errors = mode_errors
    @reset method.modes = modes
    sim.method = method

    return modes, coefficients
end


function modes_coefficients!(sim::BearingSimulation{PriorMethod})

    ω = sim.ω
    T = typeof(ω)

    bearing = sim.bearing
    method = sim.method
    
    boundarydata1 = sim.boundarydata1
    boundarydata2 = sim.boundarydata2
    
    boundarybasis1 = deepcopy(sim.boundarybasis1)
    boundarybasis2 = deepcopy(sim.boundarybasis2)

    modes = method.modal_method.modes

    ## The forward problem
        mode_errors = zeros(T, modes |> length)
        Ms = [zeros(Complex{T},4,4) for n = modes]

        is = sortperm(modes, by = abs)

        for i in is
            M = boundarycondition_system(ω, bearing, boundarybasis1.basis[1].boundarytype, boundarybasis2.basis[1].boundarytype, modes[i])
            mode_errors[i] = cond(M) * eps(T)
            Ms[i] = M

            if mode_errors[i] > method.modal_method.tol
                @warn "The relative error for the boundary conditions was $(mode_errors[i]) for (ω,mode) = $((ω,modes[i]))"
                if method.modal_method.only_stable_modes 
                    mode_errors[i] = zero(T)
                    break 
                end
            end
        end

        # just in case the loop terminated early due to keeping only_stable_modes
        i = findfirst(mode_errors[is] .== zero(T))
        if !isnothing(i)
            if i == 1
                error("there are no stable modes, and therefore no solution found")
            end

            # best to keep the indices in the same order they were given and not in the sorted order
            inds = sort(is[1:i-1])

            Ms = Ms[inds]
            mode_errors = mode_errors[inds]
            modes = modes[inds]

            # need to change all the boundarybasis modes if they exist
            basis1 = map(boundarybasis1.basis) do b
                if !isempty(b.coefficients)
                   @reset b.coefficients = b.coefficients[inds,:]
                end
                b
            end;
            boundarybasis1 = BoundaryBasis(basis1);
            
            basis2 = map(boundarybasis2.basis) do b
                if !isempty(b.coefficients)
                   @reset b.coefficients = b.coefficients[inds,:]
                end
                b
            end;
            boundarybasis2 = BoundaryBasis(basis2);
        end
        
        @reset method.modal_method.mode_errors = mode_errors
        @reset method.modal_method.modes = modes

        M_forward = BlockDiagonal(Ms)

    # calculate the prior matrix and bias vector
        prior_matrix, bias_vector = prior_and_bias(modes, boundarydata1, boundarydata2, boundarybasis1, boundarybasis2) 

          
    # Define the linear prior matrix B and prior c
        # δ = method.regularisation_parameter
        # x = [A; sqrt(δ) * I] \ [b; zeros(size(A)[2])]

        # there is not point in regularising as only well posed Ms are used
        B = M_forward \ prior_matrix
        c = M_forward \ bias_vector

        # a = vcat(wave.potentials[1].coefficients,wave.potentials[2].coefficients)[:]
        # a = B*x + c
        # M_forward * a = prior_matrix * x + bias_vector
        # norm(B * [1.0,1.0] + c - a) / norm(a)

        # prior_matrix * x 
        # M_forward * a
        # 
        # x = prior_matrix \ (M_forward * a)
        # norm(prior_matrix * x - M_forward * a) / norm(prior_matrix * x)
        # norm(prior_matrix * [1.0,1.0] - M_forward * a) / norm(prior_matrix * x)

## Calculate the inverse problem 

    # block boundary field data y_inv
    y1 = sim.boundarydata1.fields
    y2 = sim.boundarydata2.fields

    # the maximum number of measurements for any boundary
    M1 = size(sim.boundarydata1.fields,1)
    M2 = size(sim.boundarydata2.fields,1)
    
    M = max(M1,M2)

    # pad boundary data with zeros
    y1 = zeros(Complex{T},M,2)
    y2 = zeros(Complex{T},M,2)

    y1[1:M1,:] = sim.boundarydata1.fields
    y2[1:M2,:] = sim.boundarydata2.fields

    y_inv = hcat(y1, y2) |> transpose
    y_inv = y_inv[:]

    # Use the prior to solve the inverse problem
    Mns = [
        boundarycondition_system(ω, bearing, sim.boundarydata1.boundarytype, sim.boundarydata2.boundarytype, n) 
    for n in modes]
       
    M_inverse = BlockDiagonal(Mns)

    # M0_inverse = vcat([sum(M_inverse[(n+1):4:(end),:],dims = 1) for n = 0:3]...);
    # M0_inverse * a
    # M0_inverse * a - sum(Mns[i] * as[:,i] for i in eachindex(Mns))
    # abs.(M0_inverse * a - ys) 
    # abs.(sum(Mns[i] * as[:,i] for i in eachindex(Mns)) - ys) 
    # norm(sum(Mns[i] * as[:,i] for i in eachindex(Mns)) - ys) / norm(ys)

## Calculate the block matrix E where E * M_inverse * a is how the potentials contribute to the fields of the boundary conditions

    # pad angles and characteristic indicators with zeros
    θ1 = zeros(T,M); θ2 = zeros(T,M);
    χ1 = zeros(T,M); χ2 = zeros(T,M);

    θ1[1:M1] = sim.boundarydata1.θs;
    θ2[1:M2] = sim.boundarydata2.θs;

    χ1[1:M1] .= one(T)
    χ2[1:M2] .= one(T)

    Id = Matrix(one(T)I, 2, 2)
    Z = zeros(T,2,2)

    Es = [
        [Id*χ1[m]*exp(im*n*θ1[m]) Z; 
        Z Id*χ2[m]*exp(im*n*θ2[m])]
    for m = 1:M, n in modes];

    EM_inverse = Matrix(mortar(Es)) * M_inverse;

## Solve with the prior method

    # EM_inverse * (B * x + c) =  y_inv 

    # we should probably consider a regulariser here
    x = (EM_inverse * B) \ (y_inv - EM_inverse * c)

    a = B*x + c

    boundary_error = norm(EM_inverse*a - y_inv) / norm(y_inv)
    condition_number = cond(EM_inverse * B)

    @reset method.boundary_error = boundary_error
    @reset method.condition_number = condition_number
    sim.method = method

    # if relative_error > sim.method.tol
#         @warn "The relative error for the boundary conditions was $(relative_error) for (ω,basis_order) = $((ω,basis_order))"
    # end

    coefficients = a
    coefficients = reshape(coefficients,(4,:)) |> transpose |> collect

    return modes, coefficients
end

function prior_and_bias(modes::AbstractVector{Int}, boundarydata1::BoundaryData{BC1,T}, boundarydata2::BoundaryData{BC2,T}, boundarybasis1::BoundaryBasis, boundarybasis2::BoundaryBasis) where {BC1, BC2, T}

# Create the Prior matrix P 
    N1 = length(boundarybasis1.basis)
    F1s = [b.coefficients for b in boundarybasis1.basis]

    F1 = vcat((transpose.(F1s))...)
    P1s = [reshape(F1[:,m],2,N1) for m = 1:size(F1,2)]
    P1 = convert(Array{Complex{T}},vcat([ [P; 0.0*P]  for P in P1s]...))

    N2 = length(boundarybasis2.basis)
    F2s = [b.coefficients for b in boundarybasis2.basis]

    F2 = vcat((transpose.(F2s))...)
    P2s = [reshape(F2[:,m],2,N2) for m = 1:size(F2,2)]
    P2 = convert(Array{Complex{T}}, vcat([ [0.0*P; P]  for P in P2s]...))

    # reshape just in case matrix is empty so hcat can work seemlesly
    P1 = reshape(P1, 4 * length(modes), :)
    P2 = reshape(P2, 4 * length(modes), :)

    P = hcat(P1,P2)
    
# Create the Bias vector b
    # If one of the boundary basis is empty we need more boundarydata to resolve the problem. We check to see if the user has passed this boundarydata
    # bool_basis = isempty.([P1,P2])
    boundarydatas = [boundarydata1, boundarydata2]
    boundarydata_types = [b.boundarytype for b in boundarydatas]

    zs = zeros(Complex{T}, length(modes))

    b = if isempty(P1)
        i = findfirst(bt -> boundarybasis1.basis[1].boundarytype == bt, boundarydata_types)
        if isnothing(i) 
            @error "The basis boundarybasis1 was empty and there was no boundary data provided that matched the type of boundary condition given in boundarybasis1.basis[1].boundarytype. This means it is impossible to solve the forward problem."
        end

        if setdiff(modes,boundarydatas[i].modes) |> isempty
            inds = [findfirst(m .== boundarydatas[i].modes) for m in modes]
            @reset boundarydatas[i].modes = modes
            @reset boundarydatas[i].coefficients = boundarydatas[i].coefficients[inds,:]

        elseif size(boundarydatas[i].fields,1) > length(modes)
            boundarydatas[i] = fields_to_fouriermodes(boundarydatas[i],modes)
        
        else error("Not enough data in BoundaryData $i to recover the chosen modes = $modes")
        end    

        hcat(
            boundarydatas[i].coefficients[:,1], 
            boundarydatas[i].coefficients[:,2], 
            zs, 
            zs
        ) |> transpose
    else
    
        hcat(zs,zs,zs,zs) |> transpose     
    end    

    b = if isempty(P2)
        i = findfirst(bt -> boundarybasis2.basis[1].boundarytype == bt, boundarydata_types)
        if isnothing(i) 
            @error "The basis boundarybasis2 was empty and there was no boundary data provided that matched the type of boundary condition given in boundarybasis2.basis[1].boundarytype. This means it is impossible to solve the forward problem."
        end

        if setdiff(modes,boundarydatas[i].modes) |> isempty
            inds = [findfirst(m .== boundarydatas[i].modes) for m in modes]
            @reset boundarydatas[i].modes = modes
            @reset boundarydatas[i].coefficients = boundarydatas[i].coefficients[inds,:]

        elseif size(boundarydatas[i].fields,1) > length(modes)
            boundarydatas[i] = fields_to_fouriermodes(boundarydatas[i],modes)
        
        else error("Not enough data in BoundaryData $i to recover the chosen modes = $modes")
        end
    
        transpose(b) + hcat(
            zs, zs, 
            boundarydatas[i].coefficients[:,1], 
            boundarydatas[i].coefficients[:,2]
        ) |> transpose
    else b    
    end
    
    # flatten to solve the matrix system
    bias_vector = b[:]

    return P, bias_vector
end    