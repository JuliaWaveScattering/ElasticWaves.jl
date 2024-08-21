# Here we test the straight forward inverse problem for steel with a defect. See the docs/examples/steel-defect.jl for more details and plots. The example also adds errors, where here we are just checking for correctness without errors

@testset "Steel bearing defect" begin

medium = Elastic(2; ρ = 7000.0, cp = 5000.0 - 0.0im, cs = 3500.0 - 0.0im)

# need a higher rotation speed Ω for the forward and inverse problem to be well posed. In practice, we need only the inverse problem to be well posed.
Ω = 2pi * 120 / 60 # large wind turbines rotate at about 15 rpm
Z = 15 

inner_radius = 1.0
bearing = RollerBearing(medium = medium, 
    inner_radius = inner_radius, outer_radius = 1.1, 
    angular_speed = Ω,  
    rollers_inside = true,
    number_of_rollers = Z,
    roller_radius = inner_radius / (typeof(inner_radius)(1.0) + 2.5*Z / (2pi))
)

frequency_order = 5
ms = 1:frequency_order
ωms = natural_frequencies(bearing, frequency_order) |> collect

# Types of boundary conditions for the forward and inverse problem
bc1_forward = TractionBoundary(inner=true)
bc2_forward = TractionBoundary(outer=true)

bc1_inverse = DisplacementBoundary(inner=true)
bc2_inverse = TractionBoundary(inner=true)

modes_to_measure = 0:12
max_condition_number = 1e7
tol = max_condition_number * eps(Float64)

# start with stribeck, then add defects by multiplying
segment_number = 10
θs = LinRange(-pi,pi, segment_number * bearing.number_of_rollers)

q0 = 6.0; ε = 0.5;
qs = q0 .* (1 .+ 0im .- 1/(2ε) .* (1 .- cos.(θs))) .^ (10.0 / 9.0)
qs = real.(qs)
qs = [ ((q < 0) ? 0.0 : q) for q in qs]

θos = [-1.0];
# θos = [-1.0, -1.2];
is = [findmin(abs.(θo .- θs))[2] for θo in θos]
qs[is] .= 0.0

loading_profile = BoundaryData(bc1_forward, 
    θs = θs, 
    fields = hcat(qs .+ 0.0im,0.0 .* qs)
)
loading_profile = fields_to_fouriermodes(loading_profile)

## Below we solve the inverse problem

numberofsensors = 2 * length(modes_to_measure) - 1
θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]
modes_inv = -(length(modes_to_measure) - 1):(length(modes_to_measure) - 1) |> collect

bd2_for = BoundaryData(bc2_forward, 
    θs = θs, 
    fields = [zeros(Complex{Float64},  θs |> length) zeros(Complex{Float64}, θs |> length)]
)


bd2_inverse = bd2_for # a little bit of an inverse crime
modal_method = ModalMethod(tol = tol, only_stable_modes = true)

loading_datas = map(ωms) do ω

    bd1_for = BoundaryData(ω, bearing, loading_profile);
    forward_sim = BearingSimulation(ω, modal_method, bearing, bd1_for, bd2_for);

    wave = ElasticWave(forward_sim);
    # wave.method.modes
    # maximum(wave.method.mode_errors .|> abs)

    # create the data from evaluating the forward problem.
    # In practice, we would smooth the time signal to remove the higher fourier modes, which have a large amplitude, as they are due to the smooth loading. To immitate this here, we just use a fine mesh θs, and extract the number of fourier modes which our numberofsensors can measure.
    bd1_inverse = BoundaryData(bc1_inverse, bearing.inner_radius, θs, wave)
    bd1_inverse = fields_to_fouriermodes(bd1_inverse)

    # bd1 = BoundaryData(bc1_inverse, bearing.inner_radius, wave)
    # modes = intersect(bd1.modes, bd1_inverse.modes)

    # maximum(abs.(select_modes(bd1,modes).coefficients - select_modes(bd1_inverse,modes).coefficients))
    
    # norm(select_modes(bd1,modes).coefficients - select_modes(bd1_inverse,modes).coefficients) / norm(select_modes(bd1,modes).coefficients)

    # bd1 = field_modes(wave, bearing.inner_radius, bc1_inverse.fieldtype)
    
    # check forward problem worked well
    # coes = field_modes(wave, bearing.inner_radius, bc1_forward.fieldtype)
    # bd1 = select_modes(bd1_for, wave.method.modes);
    # maximum(abs.(coes - bd1.coefficients))
    # norm(coes - bd1.coefficients) / norm(coes)

    inverse_modal_method = if maximum(wave.method.modes) < maximum(modes_inv)
        ModalMethod(modes = wave.method.modes, tol = tol, only_stable_modes = true)
    else    
        ModalMethod(modes = modes_inv, tol = tol, only_stable_modes = true)
    end

    inverse_sim = BearingSimulation(ω, inverse_modal_method, bearing, bd1_inverse, bd2_inverse);

    inverse_wave = ElasticWave(inverse_sim);
    # inverse_wave.method.modes |> show
    # inverse_wave.method.mode_errors |> maximum

    predict_traction = BoundaryData(bc1_forward,  bearing.inner_radius, inverse_wave)
    # println("For ω: ", ω)
    # modes = intersect(predict_traction.modes, bd1_for.modes)
    # coes = select_modes(predict_traction,modes).coefficients
    # error = norm(coes - select_modes(bd1_for, modes).coefficients) / norm(coes)
    # println("Traction boundary error: ", error)

    ld = LoadingBoundaryData(ω,bearing,predict_traction)
    # modes = intersect(ld.modes, loading_profile.modes)
    # maximum(abs.(select_modes(loading_profile, modes).coefficients - select_modes(ld, modes).coefficients))

    modes = [-ld.modes; ld.modes]
    coes = [conj.(ld.coefficients); ld.coefficients]

    ld = BoundaryData(ld.boundarytype, 
        modes = modes, 
        coefficients = coes
    )

    # ld = LoadingBoundaryData(ω,bearing,predict_traction)
    # modes = intersect(ld.modes, loading_profile.modes)
    # maximum(abs.(select_modes(loading_profile, modes).coefficients - select_modes(ld, modes).coefficients))

    # ld
end

modes = union(vcat([l.modes for l in loading_datas]...))
coes = Array{Complex{Float64}, 2}(undef, length(modes), 2)

for data in loading_datas
    inds = [findfirst(m .== modes) for m in data.modes]
    coes[inds,:] = data.coefficients
end

loading_predict = BoundaryData(loading_datas[1].boundarytype, 
    modes = modes, 
    coefficients = coes
)

loading_predict = loading_datas[3]

modes = intersect(loading_predict.modes, loading_profile.modes)
loading_true = select_modes(loading_profile,modes)
loading_predict = select_modes(loading_predict,modes)

@test norm(loading_true.coefficients - loading_predict.coefficients) / norm(loading_true.coefficients) < 1e-8
@test maximum(abs.(loading_true.coefficients - loading_predict.coefficients)) < 1e-9

# loading_predict = fouriermodes_to_fields(loading_predict, θs)
# plot(loading_predict.θs, loading_predict.fields[:,1] .|> abs)
# plot!(θs,qs, linestyle = :dash)

end