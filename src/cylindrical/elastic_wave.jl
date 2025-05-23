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

function source_boundarycondition_mode(ω::AbstractFloat, bc::BoundaryCondition, bearing::RollerBearing, mode::Int)
    r = (bc.inner == true) ? bearing.inner_radius : bearing.outer_radius
    inner = bc.inner
    outer = bc.outer

    select_elements = [iszero(outer) iszero(inner); iszero(outer) iszero(inner)] 
    pressure_modes = select_elements .* pressure_field_mode(ω, r, bearing.medium, mode, bc.fieldtype) 

    shear_modes = select_elements .* shear_field_mode(ω, r, bearing.medium, mode, bc.fieldtype)

    return hcat(
        pressure_modes,
        shear_modes
    )
end

function source_boundarycondition_system(ω::AbstractFloat, bearing::RollerBearing, bc1::BoundaryCondition, bc2::BoundaryCondition, mode::Int)

    if bc1 == bc2
        error("Both boundary conditions are the same. It is not possible to solve the system with just one type of boundary data/condition.")
    end     

    first_boundary = source_boundarycondition_mode(ω, bc1, bearing, mode)
    second_boundary = source_boundarycondition_mode(ω, bc2, bearing, mode)

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

    φ = HelmholtzPotential(bearing.medium.cp, kP, pressure_coefficients, modes)
    ψ = HelmholtzPotential(bearing.medium.cs, kS, shear_coefficients, modes)

    return ElasticWave(ω, bearing.medium, [φ, ψ], method)
end

function ElasticWaveVector(sims::Vector{B}) where {B <: BearingSimulation{ConstantRollerSpeedMethod}}

    bearings = [s.bearing for s in sims];
    non_dims = [s.nondimensionalise for s in sims];

    if sum(non_dims) != length(sims)
        error("all BearingSimulation need to be either all non-dimensiolised or not")
    end    
    
    nondimensionalise = (sum(non_dims) == length(sims)) ? true : false

    simcopys = map(sims) do sim 
        if nondimensionalise
            nondimensionalise!(deepcopy(sim))
        else deepcopy(sim)    
        end
    end

    modes_vec, coefficients_vec = modes_coefficients!(simcopys);

    waves = map(eachindex(sims)) do i
    
        ω = simcopys[i].ω;
        method = simcopys[i].method

        kP = ω / bearings[i].medium.cp;
        kS = ω / bearings[i].medium.cs;
        ρλ2μ = bearings[i].medium.ρ * bearings[i].medium.cp^2

        coefficients = if simcopys[i].nondimensionalise
            @reset method.loading_coefficients = method.loading_coefficients .* ρλ2μ
            coefficients_vec[i] ./ kP^2

        else coefficients_vec[i]
        end

        pressure_coefficients = coefficients[:,1:2] |> transpose
        shear_coefficients = coefficients[:,3:4] |> transpose

        φ = HelmholtzPotential(bearings[i].medium.cp, kP, pressure_coefficients, modes_vec[i])
        ψ = HelmholtzPotential(bearings[i].medium.cs, kS, shear_coefficients, modes_vec[i])

        ElasticWave(ω, bearings[i].medium, [φ, ψ], method)
    end

    return waves
end

function modes_coefficients!(sim::BearingSimulation{ModalMethod})

    ω = sim.ω
    bearing = sim.bearing
    method = sim.method

    T = typeof(ω)

    modes = sim.method.modes

    coefficients = [zeros(Complex{T}, 4) for n = modes]

    mode_errors = zeros(T, length(modes))

    is = sortperm_modes(modes)

    for i in is
        A = boundarycondition_system(ω, bearing, sim.boundarydata1.boundarytype, sim.boundarydata2.boundarytype, modes[i])
        b = [
            sim.boundarydata1.coefficients[i,:];
            sim.boundarydata2.coefficients[i,:]
        ]

        SM = diagm([T(4) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
        A = A * SM

        ## NOTE: we have commented out all regularisers. The system is too small, only 4 x4, meaning that we loss 50% - 75% accuracy when regularising, which is too much. Better to not use that mode.
         
        ## solve A*x = b with svd
            # U, S, V = svd(A);
            
            ## keep only singular values that are above the threshold
            # δ = method.regularisation_parameter
            # inds = findall(S .> S[1] * δ);
            # S = S[inds]; U = U[:,inds]; V = V[:,inds];

            # Σ = Diagonal(one(T) ./ S);

            ## now A ≈ U * Diagonal(S[inds]) * V';       
            ## we construct the solution x by projecting x onto the V subspace
            # x = V * Σ * (U') * b;
            # S[end] > S[1] * eps(T) / method.tol

            # error_x = S[1] * eps(T) / S[end] 

        ## tikinov regulariser solution commented below
            # δ = method.regularisation_parameter
            # bigA = [A; sqrt(δ) * I];
            # x = bigA \ [b; zeros(size(A)[2])]

        x = A \ b    
        
        error_x = cond(A) * eps(T)
        
        error = if norm(b) > 0
            norm(A*x - b) / norm(b);
        else 
            norm(A*x - b);
        end
        
        error = max(error_x, error)

        coefficients[i] = SM * x
        mode_errors[i] = error

        if error > method.tol || isnan(error)
            @warn("The relative error for the boundary conditions was $(error) for (ω,mode) = $((ω,modes[i]))")
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
            @warn("there are no stable modes. Will return an empty solution.")
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


function modes_coefficients!(sim::BearingSimulation{P}) where P <: AbstractPriorMethod

    EM_inverse, B_for, d_for, y_inv = prior_and_bias_inverse(sim)

    ## Solve with the prior method

    # EM_inverse * (B_for * x + c_for) =  y_inv 

    # we should probably consider a regulariser here
    x = (EM_inverse * B_for) \ (y_inv - EM_inverse * d_for)

    a = B_for * x + d_for

    boundary_error = norm(EM_inverse*a - y_inv) / norm(y_inv)
    condition_number = cond(EM_inverse * B_for)

    method = sim.method

    @reset method.boundary_error = boundary_error
    @reset method.condition_number = condition_number

    sim.method = method

    # if relative_error > sim.method.tol
    #     @warn "The relative error for the boundary conditions was $(relative_error) for (ω,basis_order) = $((ω,basis_order))"
    # end

    coefficients = a
    coefficients = reshape(coefficients,(4,:)) |> transpose |> collect

    return sim.method.modal_method.modes, coefficients
end

function modes_coefficients!(sims::Vector{B}) where B <: BearingSimulation{ConstantRollerSpeedMethod}

    data = prior_and_bias_inverse.(sims);

    Bs = [d[2] for d in data];
    ds = [d[3] for d in data];

    EBs = [d[1] * d[2] for d in data];
    Eds = [d[1] * d[3] for d in data];
    y_invs = [d[4] for d in data];

    BB = vcat(EBs...);
    DD = vcat(Eds...);
    YY = vcat(y_invs...);

    # condition matrix
    SM = diagm([4.0 / sum(abs.(BB[:,j])) for j in 1:size(BB,2)])
    BBSM = BB * SM

    # x = BB \ (YY - DD);

    # Tikinov
    δ = sqrt(sims[1].method.tol * norm(YY - DD) /  maximum(size(BB)))
    bigA = [BBSM; sqrt(δ) * I];
    x = bigA \ [YY - DD; zeros(size(BB)[2])]
    
    boundary_error = norm(BBSM*x - (YY - DD)) / norm(YY - DD)
    condition_number = cond(BB)

    x = SM * x

    for sim in sims 
        method = sim.method
        @reset method.loading_coefficients = x
        @reset method.boundary_error = boundary_error
        @reset method.condition_number = condition_number
        sim.method = method
    end

    if size(BB,1) < size(BB,2)
        @error "For the constant speed roller method, we expected the system to recover the Fourier coefficients of the loading profile to be overdetermined, but it is not. Consider using more frequencies."
    end     

    coefficients_vec = map(eachindex(sims)) do i
        coes = Bs[i] * x + ds[i]
        reshape(coes,(4,:)) |> transpose |> collect
    end    
    modes_vec = [sims[i].method.modal_method.modes for i in eachindex(sims)]

    return modes_vec, coefficients_vec
end

function prior_and_bias_inverse(sim::BearingSimulation{P}) where P <: AbstractPriorMethod

    ω = sim.ω
    T = typeof(ω)

    bearing = sim.bearing
    method = sim.method
    
    boundarydata1 = sim.boundarydata1
    boundarydata2 = sim.boundarydata2
    
    boundarybasis1 = deepcopy(sim.boundarybasis1)
    boundarybasis2 = deepcopy(sim.boundarybasis2)

    method_modes = method.modal_method.modes

    ## The forward problem
        mode_errors = zeros(T, method_modes |> length)
        MSs = [zeros(Complex{T},4,4) for n = method_modes]
        Ss = [zeros(Complex{T},4,4) for n = method_modes]

        # this sortperm_modes below is now depricated, as all modes arrive at this point already sorted in this order. Will remove 
        is = sortperm_modes(method_modes)

        for i in is
            M = boundarycondition_system(ω, bearing, boundarybasis1.basis[1].boundarytype, boundarybasis2.basis[1].boundarytype, method_modes[i])

            S = diagm([T(4) / sum(abs.(M[:,j])) for j in 1:size(M,2)])
            MS = M * S

            mode_errors[i] = cond(MS) * eps(T)
            MSs[i] = MS
            Ss[i] = S

            if mode_errors[i] > method.modal_method.tol
                @warn "The relative error for the boundary conditions was $(mode_errors[i]) for (ω,mode) = $((ω,method_modes[i]))"
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

            MSs = MSs[inds]
            Ss = Ss[inds]
            mode_errors = mode_errors[inds]
            method_modes = method_modes[inds]

            # need to change all the boundarybasis modes if they exist
            basis1 = map(boundarybasis1.basis) do b
                if !isempty(b.coefficients)
                   @reset b.coefficients = b.coefficients[inds,:]
                   @reset b.modes = b.modes[inds]
                end
                b
            end;
            boundarybasis1 = BoundaryBasis(basis1);
            
            basis2 = map(boundarybasis2.basis) do b
                if !isempty(b.coefficients)
                   @reset b.coefficients = b.coefficients[inds,:]
                   @reset b.modes = b.modes[inds]
                end
                b
            end;
            boundarybasis2 = BoundaryBasis(basis2);
        end
        
        @reset method.modal_method.mode_errors = mode_errors
        @reset method.modal_method.modes = method_modes

        sim.method = method

        M_forward = BlockDiagonal(MSs)
        S_forward = BlockDiagonal(Ss)

    # calculate the prior matrix and bias vector
        prior_matrix, bias_vector = prior_and_bias_forward(method_modes, boundarydata1, boundarydata2, boundarybasis1, boundarybasis2) 
          
    # Define the linear prior matrix B and prior c
        # δ = method.regularisation_parameter
        # x = [A; sqrt(δ) * I] \ [b; zeros(size(A)[2])]

        # there is no point in regularising as only well posed Ms are used
        B_for = S_forward * (M_forward \ prior_matrix)
        c_for = S_forward * (M_forward \ bias_vector)

        # a = vcat(wave.potentials[1].coefficients,wave.potentials[2].coefficients)[:]
        # a = B_for*x + c_for
        # inv(S_forward) * M_forward * a = prior_matrix * x + bias_vector
        # norm(B_for * [1.0,1.0] + c - a) / norm(a)

        # prior_matrix * x 
        # M_forward * a
        # 
        # x = prior_matrix \ (inv(S_forward) * M_forward * a)
        # norm(prior_matrix * x - M_forward * a) / norm(prior_matrix * x)
        # norm(prior_matrix * [1.0,1.0] - M_forward * a) / norm(prior_matrix * x)

## Calculate the inverse problem 

    # block boundary field data y_inv
    y1 = sim.boundarydata1.fields
    y2 = sim.boundarydata2.fields

    y_inv = vcat(transpose(y1)[:], transpose(y2)[:])

    # Use the prior to solve the inverse problem
    Mns = [
        boundarycondition_system(ω, bearing, sim.boundarydata1.boundarytype, sim.boundarydata2.boundarytype, n) 
    for n in method_modes];
       
    # M_inverse = BlockDiagonal(Mns)

    # M0_inverse = vcat([sum(M_inverse[(n+1):4:(end),:],dims = 1) for n = 0:3]...);
    # M0_inverse * a
    # M0_inverse * a - sum(Mns[i] * as[:,i] for i in eachindex(Mns))
    # abs.(M0_inverse * a - ys) 
    # abs.(sum(Mns[i] * as[:,i] for i in eachindex(Mns)) - ys) 
    # norm(sum(Mns[i] * as[:,i] for i in eachindex(Mns)) - ys) / norm(ys)

## Calculate the block matrix E where E * M_inverse * a is how the potentials contribute to the fields of the boundary conditions

    θ1 = sim.boundarydata1.θs;
    θ2 = sim.boundarydata2.θs;

    E1s = [
        exp(im * method_modes[j] * θ) .* Mns[j][1:2,:]
    for θ in θ1, j in eachindex(method_modes)];
    EE1 = Matrix(mortar(Matrix{Matrix{Complex{T}}}(E1s)));# This is used to deal with an empty matrix
    EE1 = reshape(EE1, :, length(method_modes) * size(Mns[1],2));

    E2s = [
        exp(im * method_modes[j] * θ) .* Mns[j][3:4,:]
    for θ in θ2, j in eachindex(method_modes)];
    EE2 = Matrix(mortar(Matrix{Matrix{Complex{T}}}(E2s)));# This is used to deal with an empty matrix
    EE2 = reshape(EE2, :, length(method_modes) * size(Mns[1],2));

    EM_inverse = vcat(EE1,EE2);

    return EM_inverse, B_for, c_for, y_inv
end

function prior_and_bias_forward(modes::AbstractVector{Int}, boundarydata1::BoundaryData{BC1,T}, boundarydata2::BoundaryData{BC2,T}, boundarybasis1::BoundaryBasis, boundarybasis2::BoundaryBasis) where {BC1, BC2, T}

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