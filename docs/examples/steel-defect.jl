
# Here we test 

using ElasticWaves
using Test, Statistics, LinearAlgebra, MultipleScattering
using Plots

medium = Elastic(2; ρ = 7000.0, cp = 5000.0 - 0.0im, cs = 3500.0 - 0.0im)

# need a higher rotation speed Ω for the forward and inverse problem to be well posed. In practice, we need only the inverse problem to be well posed.
Ω = 2pi * 120 / 60 # large wind turbines rotate at about 15 rpm
# Ω = 2pi * 1200 / 60 # large wind turbines rotate at about 15 rpm
Z = 15

# c =  

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
modes_to_measure = 0:20
max_condition_number = 1e7
max_condition_number = 1e5
max_condition_number = 2e8
tol = max_condition_number * eps(Float64)

ωs = ωms
modes_vec = map(ωms) do ω
    conds = map(modes_to_measure) do n
        nondim_bearing = nondimensionalise(bearing, ω)
        A = boundarycondition_system(ω, nondim_bearing, bc1_forward, bc2_forward, n) 
        # A = boundarycondition_system(ω, nondim_bearing, bc1_inverse, bc2_inverse, n) 
        SM = diagm([(4.0) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
        A = A * SM
        cond(A)
    end
    inds = findall(conds .< max_condition_number)
    modes_to_measure[inds]
end

loading_modes_vec = map(eachindex(ms)) do i
    modes_vec[i] .- ms[i] * Z
end

scatter_ms_vec = map(eachindex(ms)) do i
    [ms[i] for l in loading_modes_vec[i]]
end
scatter(vcat(scatter_ms_vec...),abs.(vcat(loading_modes_vec...)))
scatter!(ylims = (0, 50), xlab = "ωm frequency", ylab = "loading modes")

measurable_loading_modes = union(vcat(loading_modes_vec...));
all_measurable_modes = [-reverse(measurable_loading_modes); measurable_loading_modes]

# start with stribeck, then add defects by multiplying
segment_number = 10
θs = LinRange(-pi,pi, segment_number * bearing.number_of_rollers)[1:end-1]

q0 = 6.0; ε = 0.5;
qs = q0 .* (1 .+ 0im .- 1/(2ε) .* (1 .- cos.(θs))) .^ (10.0 / 9.0)
qs = real.(qs)
qs = [ ((q < 0) ? 0.0 : q) for q in qs]

θos = [-1.0];
θos = [-1.0, 0.1];
is = [findmin(abs.(θo .- θs))[2] for θo in θos]

qs_defect = θs .* 0.0 |> collect;
qs_defect[is] = - qs[is]

qs = qs + qs_defect 
plot(θs,qs)

loading_profile = BoundaryData(bc1_forward, 
    θs = θs, 
    fields = hcat(qs .+ 0.0im,0.0 .* qs)
)
loading_profile = fields_to_fouriermodes(loading_profile)
bd1_for = BoundaryData(ωms[3], bearing, loading_profile);
# scatter(loading_profile.modes, loading_profile.coefficients[:,1] .|> abs, ylims = (0.,0.1))
scatter(bd1_for.modes, bd1_for.coefficients[:,1] .|> abs, ylims = (0.,0.2))


non_measurable_modes = setdiff(loading_profile.modes, all_measurable_modes)

## Below we assume the inverse problem will work where it is well conditioned and just check what the expected recovery looks like

    zero_inds = [
        findfirst(m .== loading_profile.modes) 
    for m in non_measurable_modes]

    loading_profile2 = deepcopy(loading_profile)
    loading_profile2.coefficients[zero_inds,:] .= 0.0 + 0.0im
    loading_profile2 = fouriermodes_to_fields(loading_profile2)

    inds = sortperm(loading_profile.modes)
    loading_profile.modes[inds]

    fun = real
    plot(loading_profile.modes[inds],loading_profile.coefficients[inds,1] .|> fun)
    scatter!(loading_profile.modes[inds],loading_profile2.coefficients[inds,1] .|> fun)
    plot!(ylims = (-0.04,0.04), xlims = (0,75))

    fun = abs
    plot(loading_profile.θs, loading_profile.fields .|> fun)
    plot!(loading_profile2.θs, loading_profile2.fields .|> fun, linestyle = :dash)

## Below we solve the inverse problem
bd2_for = BoundaryData(bc2_forward, 
    θs = θs, 
    fields = [zeros(Complex{Float64},  θs |> length) zeros(Complex{Float64}, θs |> length)]
)

bd2_inverse = bd2_for # a little bit of an inverse crime
modal_method = ModalMethod(tol = tol, only_stable_modes = true)

function predict_boundary(rel_error = 0.005, modes_to_measure = modes_to_measure)

    numberofsensors = 2 * length(modes_to_measure) - 1
    θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]
    modes_inv = -(length(modes_to_measure) - 1):(length(modes_to_measure) - 1) |> collect

    loading_datas = map(ωms) do ω

        bd1_for = BoundaryData(ω, bearing, loading_profile);
        forward_sim = BearingSimulation(ω, modal_method, bearing, bd1_for, bd2_for);
        wave = ElasticWave(forward_sim);

        bd1_inverse = BoundaryData(bc1_inverse, bearing.inner_radius, θs, wave)

        # add rel_error
        amp = mean(abs.(bd1_inverse.fields))
        error = rel_error .* amp .* (rand(Complex{Float64},size(bd1_inverse.fields)...) .- 0.5 .- 0.5im)

        bd1_inverse = BoundaryData(bc1_inverse, 
            θs = θs, 
            fields = bd1_inverse.fields + error
        )
        bd1_inverse = fields_to_fouriermodes(bd1_inverse)

        inverse_modal_method = if maximum(wave.method.modes) < maximum(modes_inv)
            ModalMethod(modes = wave.method.modes, tol = tol, only_stable_modes = true)
        else    
            ModalMethod(modes = modes_inv, tol = tol, only_stable_modes = true)
        end

        inverse_sim = BearingSimulation(ω, inverse_modal_method, bearing, bd1_inverse, bd2_inverse);
        inverse_wave = ElasticWave(inverse_sim);

        predict_traction = BoundaryData(bc1_forward,  bearing.inner_radius, inverse_wave)
        ld = LoadingBoundaryData(ω,bearing,predict_traction)

        modes = [-ld.modes; ld.modes]
        coes = [conj.(ld.coefficients); ld.coefficients]

        ld = BoundaryData(ld.boundarytype, 
            modes = modes, 
            coefficients = coes
        )
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

    return fouriermodes_to_fields(loading_predict, θs)
end

measured_modes = 0:15
loading_predict = predict_boundary(0.01,measured_modes)

modes = intersect(loading_predict.modes, loading_profile.modes)
loading_true = select_modes(loading_profile,modes)
loading_predict = select_modes(loading_predict,modes)

norm(loading_true.coefficients - loading_predict.coefficients) / norm(loading_true.coefficients)

loading_predict = fouriermodes_to_fields(loading_predict, θs)

plot(loading_predict.θs, loading_predict.fields[:,1] .|> real)
plot!(θs,qs_defect, linestyle = :dash)

norm(loading_predict.fields[:,1]- qs_defect) / norm(loading_predict.fields[:,1])

# Do a ribbon plot 

data = map(1:150) do i
    real.(predict_boundary(0.01,measured_modes).fields[:,1])
end

## redo the simulation 100 times each with a different error

# calculate the mean and standard deviation of y for each value of x
y_mat = hcat(data...) # simpler to make a big matrix
ys_mean = mean(y_mat, dims = 2)[:]
ys_std = std(y_mat, dims = 2)[:]

h = 400
gr(linewidth = 1.0, size = (1.6 * h, h ))

plot(θs, abs.(ys_mean); ribbon = ys_std .* 0.5, fillalpha=.4, linewidth = 2.0, lab = "predicted")
plot!(θs,qs, linestyle = :dash, lab = "exact")
plot!(xlab = "θ", ylab = "loading profile")

# savefig("docs/images/steel-bearing-defect-modes-$(measured_modes[end])-RPM-$(Int(round(60 * Ω / (2pi)))).pdf")

# use fewer modes
measured_modes = 0:6

loading_predict = predict_boundary(0.01,measured_modes)

modes = intersect(loading_predict.modes, loading_profile.modes)
loading_true = select_modes(loading_profile,modes)
loading_predict = select_modes(loading_predict,modes)

norm(loading_true.coefficients - loading_predict.coefficients) / norm(loading_true.coefficients)

loading_predict = fouriermodes_to_fields(loading_predict, θs)

plot(loading_predict.θs, loading_predict.fields[:,1] .|> real)
plot!(θs,qs_defect, linestyle = :dash)

norm(loading_predict.fields[:,1]- qs_defect) / norm(loading_predict.fields[:,1])

# Do a ribbon plot 

data = map(1:100) do i
    real.(predict_boundary(0.02,measured_modes).fields[:,1])
end

## redo the simulation 100 times each with a different error

# calculate the mean and standard deviation of y for each value of x
y_mat = hcat(data...) # simpler to make a big matrix
ys_mean = mean(y_mat, dims = 2)[:]
ys_std = std(y_mat, dims = 2)[:]

h = 250
gr(linewidth = 1.0, size = (1.6 * h, h ))

plot(θs, abs.(ys_mean); ribbon = ys_std .* 0.5, fillalpha=.4, linewidth = 2.0, lab = "predicted")
plot!(θs,qs, linestyle = :dash, lab = "exact")
plot!(xlab = "θ", ylab = "loading profile")
# plot!(xlab = "θ", ylab = "loading profile", ylims = (-2., 6.))

# savefig("docs/images/steel-bearing-defect-6-modes-error-2%.pdf")