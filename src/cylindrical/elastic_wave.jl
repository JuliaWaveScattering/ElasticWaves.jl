function boundarycondition_mode(ω::AbstractFloat, bc::BoundaryCondition, bearing::RollerBearing, basis_order::Int)
    r = (bc.inner == true) ? bearing.inner_radius : bearing.outer_radius
    return hcat(
        pressure_field_mode(ω, r, bearing.medium, basis_order, bc.fieldtype),
        shear_field_mode(ω, r, bearing.medium, basis_order, bc.fieldtype)
    )
end

function boundarycondition_system(ω::AbstractFloat, bearing::RollerBearing, bc1::BoundaryCondition, bc2::BoundaryCondition, basis_order::Int)

    if bc1 == bc2
        error("Both boundary conditions are the same. It is not possible to solve the system with just one type of boundary data/condition.")
    end     

    first_boundary = boundarycondition_mode(ω, bc1, bearing, basis_order)
    second_boundary = boundarycondition_mode(ω, bc2, bearing, basis_order)

    return vcat(first_boundary, second_boundary)
end

# replace the below with
# function ElasticWave(sim::BearingSimulation{ModalMethod})

# write another function for the priors
# function ElasticWave(sim::BearingSimulation{PriorMethod})

function ElasticWave(sim::BearingSimulation{ModalMethod})

    ω = sim.ω
    bearing = sim.bearing

    T = typeof(ω)

    kP = ω / bearing.medium.cp;
    kS = ω / bearing.medium.cs

    basis_order = sim.basis_order

    coefficients = [
        zeros(Complex{T}, 4)
    for m = -basis_order:basis_order]

    mode_errors = zeros(T, 2basis_order + 1)

    for m in -basis_order:basis_order
        A = boundarycondition_system(ω, bearing, sim.boundarydata1.boundarytype, sim.boundarydata2.boundarytype, m)
        b = [
             sim.boundarydata1.fourier_modes[m+basis_order+1,:];
             sim.boundarydata2.fourier_modes[m+basis_order+1,:]
        ]

        # solve A*x = b with tikinov regulariser
        # x = A \ b
        δ = sim.method.regularisation_parameter
        x = [A; sqrt(δ) * I] \ [b; zeros(size(A)[2])]

        relative_error = norm(A*x - b) / norm(b)
        if relative_error > sim.method.tol
            @warn "The relative error for the boundary conditions was $(relative_error) for (ω,basis_order) = $((ω,m))"
        end

        coefficients[m + basis_order + 1] = x
        mode_errors[m + basis_order + 1] = relative_error
    end

    coefficients = transpose(hcat(coefficients...)) |> collect

    if sim.method.only_stable_modes
        inds = findall(mode_errors .< sim.method.tol)

        orders = inds .- (basis_order + 1)

        # find a contiguous block with small mode_errors
        m = min(abs(orders[1]),abs(orders[end]))
        unstable_orders = setdiff(-m:m,orders)
        if !(unstable_orders |> isempty)
            m = minimum(abs.(unstable_orders)) - 1
            if m < 0 
                @error "no stable modes found for the given tolerance tol = $(sim.method.tol)"
            end 
        end

        inds = (-m:m) .+ (basis_order + 1)
        
        coefficients = coefficients[inds,:]
        mode_errors = mode_errors[inds]
    end

    pressure_coefficients = coefficients[:,1:2] |> transpose
    shear_coefficients = coefficients[:,3:4] |> transpose

    φ = HelmholtzPotential{2}(bearing.medium.cp, kP, pressure_coefficients)
    ψ = HelmholtzPotential{2}(bearing.medium.cs, kS, shear_coefficients)

    return ElasticWave(ω, bearing.medium, [φ, ψ]; mode_errors = mode_errors)
end


function ElasticWave(sim::BearingSimulation{PriorMethod})

    ω = sim.ω
    bearing = sim.bearing
    
    boundarydata1 = sim.boundarydata1
    boundarydata2 = sim.boundarydata2
    
    boundarybasis1 = sim.boundarybasis1
    boundarybasis2 = sim.boundarybasis2

    basis_order = sim.basis_order

    # check which orders are stable for the forward problem
        mode_errors = zeros(T, 2basis_order + 1)

        Ms = [
            boundarycondition_system(ω, bearing, boundarybasis1.basis[1].boundarytype, boundarybasis2.basis[1].boundarytype, n)  
        for n in -basis_order:basis_order]

        M_forward = BlockDiagonal(Ms)


    # Create the Prior matrix P 
        N1 = length(boundarybasis1.basis)
        F1s = [b.fourier_modes for b in boundarybasis1.basis]

        F1 = vcat((transpose.(F1s))...)
        P1s = [reshape(F1[:,m],2,N1) for m = 1:size(F1,2)]
        P1 = vcat([ [P; 0.0*P]  for P in P1s]...)

        N2 = length(boundarybasis2.basis)
        F2s = [b.fourier_modes for b in boundarybasis2.basis]

        F2 = vcat((transpose.(F2s))...)
        P2s = [reshape(F2[:,m],2,N2) for m = 1:size(F2,2)]
        P2 = vcat([ [0.0*P; P]  for P in P2s]...)

        # reshape just in case matrix is empty so hcat can work seemlesly
        P1 = reshape(P1, 4 * (2 * sim.basis_order + 1), :)
        P2 = reshape(P2, 4 * (2 * sim.basis_order + 1), :)
    
        P = hcat(P1,P2)
        
    # Create the Bias vector b
        # If one of the boundary basis is empty we need more boundarydata to resolve the problem. We check to see if the user has passed this boundarydata
        # bool_basis = isempty.([P1,P2])
        boundarydatas = [boundarydata1, boundarydata2]
        boundarydata_types = [b.boundarytype for b in boundarydatas]

        b = if isempty(P1)
            i = findfirst(bt -> boundarybasis1.basis[1].boundarytype == bt, boundarydata_types)
            if isnothing(i) 
                @error "The basis boundarybasis1 was empty and there was no boundary data provided that matched the type of boundary condition given in boundarybasis1.basis[1].boundarytype. This means it is impossible to solve the forward problem."
            end

            hcat(
                boundarydatas[i].fourier_modes[:,1], 
                boundarydatas[i].fourier_modes[:,2], 
                0.0 .* boundarydatas[i].fourier_modes[:,1], 
                0.0 .* boundarydatas[i].fourier_modes[:,1]
            ) |> transpose
        else
            b = hcat(
                0.0 .* boundarydatas[1].fourier_modes[:,1], 
                0.0 .* boundarydatas[1].fourier_modes[:,2], 
                0.0 .* boundarydatas[1].fourier_modes[:,1], 
                0.0 .* boundarydatas[1].fourier_modes[:,1]
            ) |>transpose     
        end    

        if isempty(P2)
            i = findfirst(bt -> boundarybasis2.basis[1].boundarytype == bt, boundarydata_types)
            if isnothing(i) 
                @error "The basis boundarybasis2 was empty and there was no boundary data provided that matched the type of boundary condition given in boundarybasis2.basis[1].boundarytype. This means it is impossible to solve the forward problem."
            end

            b = transpose(b) + hcat(0.0 .* boundarydatas[i].fourier_modes[:,1], 0.0 .* boundarydatas[i].fourier_modes[:,2], boundarydatas[i].fourier_modes[:,1], boundarydatas[i].fourier_modes[:,2]) |> transpose
        end
        
        # flatten to solve the matrix system
        bias_vector = b[:]
        
    T = typeof(ω)
    
    P = Matrix{ComplexF64}(P)
    
    # Define the linear prior matrix B and prior c
    B = M_forward \ P
    c = M_forward \ bias_vector

    # Use the prior to solve the inverse problem
    Ms_inverse = [
        boundarycondition_system(ω, bearing, sim.boundarydata1.boundarytype, sim.boundarydata2.boundarytype, n) 
    for n in -basis_order:basis_order]
       
    M_inverse = BlockDiagonal(Ms_inverse)

## Calculate the block boundary field data y_inv

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
    for m = 1:M, n in -basis_order:basis_order]

    EM_inverse = Matrix(mortar(Es)) * M_inverse

## Solve with the prior method

    x = (EM_inverse * B) \ (y_inv - EM_inverse * c)

    a = B*x + c

    relative_error = norm(EM_inverse*a - y_inv) / norm(y_inv)
    if relative_error > sim.method.tol
#         @warn "The relative error for the boundary conditions was $(relative_error) for (ω,basis_order) = $((ω,basis_order))"
    end

    coefficients = a
    mode_errors = relative_error

    coefficients = reshape(coefficients,(4,:)) |> transpose |> collect

    pressure_coefficients = coefficients[:,1:2] |> transpose
    shear_coefficients = coefficients[:,3:4] |> transpose

    φ = HelmholtzPotential{2}(bearing.medium.cp, kP, pressure_coefficients)
    ψ = HelmholtzPotential{2}(bearing.medium.cs, kS, shear_coefficients)

    return ElasticWave(ω, bearing.medium, [φ, ψ])

    # return ElasticWave(ω, bearing.medium, [φ, ψ]; mode_errors = mode_errors)
end
