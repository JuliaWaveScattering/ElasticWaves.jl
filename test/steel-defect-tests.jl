
# Here we test 

using ElasticWaves
using Test, Statistics, LinearAlgebra, MultipleScattering
using Accessors

# @testset "Steel bearing" begin

medium = Elastic(2; ρ = 2100.0, cp = 1500.0 - 0.0im, cs = 1050.0 - 0.0im)
medium = Elastic(2; ρ = 7000.0, cp = 5000.0 - 0.0im, cs = 3500.0 - 0.0im)

Ω = 2pi * 90 / 60 # large wind turbines rotate at about 15 rpm
Z = 15 

inner_radius = 1.0
bearing = RollerBearing(medium = medium, 
    inner_radius = inner_radius, outer_radius = 1.1, 
    angular_speed = Ω,  
    rollers_inside = true,
    number_of_rollers = Z,
    roller_radius = inner_radius / (typeof(inner_radius)(1.0) + 2.5*Z / (2pi))
)

using Plots

plot(bearing)

frequency_order = 5
ms = 1:frequency_order
ωms = natural_frequencies(bearing, frequency_order) |> collect

# estimate the number C in the restriction |n| < C ω for the system to be well conditioned. We choose on ω and increase n until we get an ill-conditioned system.

bc1_forward = TractionBoundary(inner=true)
bc2_forward = TractionBoundary(outer=true)

bc1_inverse = DisplacementBoundary(inner=true)
bc2_inverse = TractionBoundary(inner=true)

modes = 0:14
max_condition_number = 1e7
tol = max_condition_number * eps(Float64)

ωs = ωms
modes_vec = map(ωms) do ω
    conds = map(modes) do n
        nondim_bearing = nondimensionalise(bearing, ω)
        # A = boundarycondition_system(ω, nondim_bearing, bc1_forward, bc2_forward, n) 
        A = boundarycondition_system(ω, nondim_bearing, bc1_inverse, bc2_inverse, n) 
        SM = diagm([(4.0) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
        A = A * SM
        cond(A)
    end
    inds = findall(conds .< max_condition_number)
    modes[inds]
end

conds_vec = map(ωms) do ω
    map(modes) do n
        nondim_bearing = nondimensionalise(bearing, ω)
        # A = boundarycondition_system(ω, nondim_bearing, bc1_forward, bc2_forward, n) 
        A = boundarycondition_system(ω, nondim_bearing, bc1_inverse, bc2_inverse, n) 
        SM = diagm([(4.0) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
        A = A * SM
        cond(A)
    end
end

loading_modes_vec = map(eachindex(ms)) do i
    modes_vec[i] .- ms[i] * Z
end

min_loading_ns = [minimum(abs.(modes)) for modes in loading_modes_vec]
max_loading_ns = [maximum(abs.(modes)) for modes in loading_modes_vec]

scatter(min_loading_ns)
scatter!(max_loading_ns)
scatter!(ylims = (0,50), xlab = "ωm frequency", ylab = "n from loading")

scatter_ms_vec = map(eachindex(ms)) do i
    [ms[i] for l in loading_modes_vec[i]]
end

scatter(vcat(scatter_ms_vec...),abs.(vcat(loading_modes_vec...)))
scatter!(ylims = (0, 50), xlab = "ωm frequency", ylab = "loading modes")

measurable_loading_modes = union(vcat(loading_modes_vec...));
all_measurable_modes = [-reverse(measurable_loading_modes); measurable_loading_modes]

# Cs = ns ./ ωs
# min_loading_ns = (1 .- Cs .* Ω) .* ms .* Z 
# max_loading_ns = (1 .+ Cs .* Ω) .* ms .* Z 
# min_loading_ns = Int.(round.(min_loading_ns))
# max_loading_ns = Int.(round.(max_loading_ns))
# mode_intervals = hcat(min_loading_ns, max_loading_ns) |> transpose |> collect
# all_measurable_modes = [mode_intervals[1,i]:mode_intervals[2,i] for i in 1:size(mode_intervals,2)]
# all_measurable_modes = vcat(all_measurable_modes...)
# all_measurable_modes = union(all_measurable_modes)
# all_measurable_modes = [-reverse(all_measurable_modes); all_measurable_modes]

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
plot(θs,qs)

loading_profile = BoundaryData(bc1_forward, 
    θs = θs, 
    fields = hcat(qs .+ 0.0im,0.0 .* qs)
)
loading_profile = fields_to_fouriermodes(loading_profile)

non_measurable_modes = setdiff(loading_profile.modes, all_measurable_modes)

# non_measurable_modes = -12:12

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
plot!(loading_profile.θs, loading_profile2.fields .|> fun, linestyle = :dash)

## Generate sensor displacement data
numberofsensors = 100
forward_θs = θs
θs_inv = LinRange(0, 2pi, numberofsensors + 1)[1:end-1]

bd2_for = BoundaryData(bc2_forward, 
    θs = forward_θs, 
    fields = [zeros(Complex{Float64},  forward_θs |> length) zeros(Complex{Float64}, forward_θs |> length)]
)

# a little bit of an inverse crime
bd2_inverse = bd2_for
modal_method = ModalMethod(tol = tol, only_stable_modes = true)

loading_datas = map(ωms) do ω

    bd1_for = BoundaryData(ω, bearing, loading_profile);
    forward_sim = BearingSimulation(ω, modal_method, bearing, bd1_for, bd2_for);

    wave = ElasticWave(forward_sim);
    # wave.method.modes
    # maximum(wave.method.mode_errors .|> abs)

    # create the data from evaluating the forward problem 
    bd1_inverse = BoundaryData(bc1_inverse, bearing.inner_radius, θs_inv, wave)

    # to immitate smoothing in time, we should finely samply in θs and then smooth the data in the same way. Then sample only for the number of sensors


    # amp = 5e-4 .* mean(abs.(bd1_inverse.fields))
    # amp = 0 .* mean(abs.(bd1_inverse.fields))

    # reals = amp/2 .* (rand(size(bd1_inverse.fields)...) .- 0.5)
    # imags = amp/2 .* (rand(size(bd1_inverse.fields)...) .- 0.5) .* im

    # @reset bd1_inverse.fields = bd1_inverse.fields + reals + imags 

    inverse_modal_method = if wave.method.modes |> length < length(θs_inv)
        ModalMethod(modes = wave.method.modes,tol = tol, only_stable_modes = true)
    else    
        ModalMethod(tol = tol, only_stable_modes = true)
    end

    inverse_sim = BearingSimulation(ω, inverse_modal_method, bearing, bd1_inverse, bd2_inverse);

    inverse_wave = ElasticWave(inverse_sim);
    inverse_wave.method.modes |> show
    inverse_wave.method.mode_errors |> show

    predict_traction = BoundaryData(bc1_forward,  bearing.inner_radius, inverse_wave)

    # norm(predict_traction.coefficients - select_modes(bd1_for, predict_traction.modes).coefficients) / norm(predict_traction.coefficients)

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

loading_data = BoundaryData(loading_datas[1].boundarytype, 
    modes = modes, 
    coefficients = coes
)

loading_data = fouriermodes_to_fields(loading_data, θs)

plot(loading_data.θs, 1.85 .+ loading_data.fields[:,1] .|> real)
plot!(θs,qs, linestyle = :dash)


loading_data = loading_datas[4]
modes = intersect(loading_data.modes, loading_profile.modes)
loading_data = select_modes(loading_data, modes)

ld_true = select_modes(loading_profile,loading_data.modes)

maximum(abs.(ld_true.coefficients - loading_data.coefficients))
norm(ld_true.coefficients - loading_data.coefficients) / norm(loading_data.coefficients)