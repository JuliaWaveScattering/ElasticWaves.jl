# Here we test 

# @testset "Steel bearing" begin

# Well almost steel!
# medium = Elastic(2; ρ = 2.0, cp = 5000.0 - 0.0im, cs = 3500.0 - 0.0im)
medium = Elastic(2; ρ = 2.0, cp = 500.0 - 0.0im, cs = 350.0 - 0.0im)

Ω = 2pi * 15 / 60 # large wind turbines rotate at about 15 rpm
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

frequency_order = 6
ms = 1:frequency_order
ωms = natural_frequencies(bearing, frequency_order) |> collect

# estimate the number C in the restriction |n| < C ω for the system to be well conditioned. We choose on ω and increase n until we get an ill-conditioned system.

# For k a = 1 we need
# kp = ω / real(medium.cp)

bc1_forward = TractionBoundary(inner=true)
bc2_forward = TractionBoundary(outer=true)
ns = 0:20
# ωs = (real(medium.cp) / inner_radius) .* [0.01,0.1,1.0]

# min_condition_number = 4e-5
max_condition_number = 1e7
ωs = ωms
Cs = map(ωs) do ω
    conds = map(ns) do n
        A = boundarycondition_system(ω, bearing, bc1_forward, bc2_forward, n) 
        SM = diagm([(4.0) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
        A = A * SM

        # A |> det |> abs
        cond(A)
    end
    # plot(ns,conds, ylims = (0,0.001))

    # n = findfirst(conds .< min_condition_number) - 1
    n = findfirst(conds .> max_condition_number)
    n = isnothing(n) ? ns[end] : n 
    n / ω
end

min_loading_ns = (1 .- Cs .* Ω) .* ms .* Z 
max_loading_ns = (1 .+ Cs .* Ω) .* ms .* Z 

mode_intervals = hcat(min_loading_ns, max_loading_ns) |> transpose |> collect

plot(min_loading_ns)
plot!(max_loading_ns)

# start with stribeck, then add defects by multiplying
segment_number = 10
θs = LinRange(-pi,pi, segment_number * bearing.number_of_rollers)

q0 = 6.0; ε = 0.5;
qs = q0 .* (1 .+ 0im .- 1/(2ε) .* (1 .- cos.(θs))) .^ (10.0 / 9.0)
qs = real.(qs)
qs = [ ((q < 0) ? 0.0 : q) for q in qs]

θo = -1.0;
i = findmin(abs.(θo .- θs))[2]
qs[i] = 0.0

# θo = -1.2;
# i = findmin(abs.(θo .- θs))[2]
# qs[i] = 0.0


plot(θs,qs)

loading_profile = BoundaryData(bc1_forward, 
    θs = θs, 
    fields = hcat(qs .+ 0.0im,0.0 .* qs)
)

loading_profile = fields_to_fouriermodes(loading_profile)


all_measurable_modes = [mode_intervals[1,i]:mode_intervals[2,i] for i in 1:size(mode_intervals,2)]

all_measurable_modes = vcat(all_measurable_modes...)
all_measurable_modes = union(all_measurable_modes)

all_measurable_modes = [-reverse(all_measurable_modes); all_measurable_modes]

non_measurable_modes = setdiff(loading_profile.modes, all_measurable_modes)

# non_measurable_modes = -12:12

zero_inds = [
    findfirst(m .== loading_profile.modes) 
for m in non_measurable_modes]

measured_inds = [
    findfirst(m .== loading_profile.modes) 
for m in all_measurable_modes]

inds = findall(measured_inds .== nothing)
deleteat!(measured_inds, inds)
measured_inds = Int.(measured_inds)


loading_profile2 = deepcopy(loading_profile)
loading_profile2.coefficients[zero_inds,:] .= 0.0 + 0.0im

inds = sortperm(loading_profile.modes)
loading_profile.modes[inds]

fun = real
plot(loading_profile.modes[inds],loading_profile.coefficients[inds,1] .|> fun)
plot!(ylims = (-0.04,0.04), xlims = (10,70))
scatter!(loading_profile.modes[inds],loading_profile2.coefficients[inds,1] .|> fun)
plot!(ylims = (-0.04,0.04), xlims = (10,70))

using DataInterpolations
 
coefficients = loading_profile2.coefficients[measured_inds,:] 
modes = loading_profile2.modes[measured_inds] 

inds = sortperm(modes)

coes = CubicSpline(coefficients[inds,1], modes[inds])

abs_coes = CubicSpline(abs.(coefficients[inds,1]), modes[inds])
angle_coes = CubicSpline(angle.(coefficients[inds,1]), modes[inds])
# coes = QuadraticSpline(coefficients[inds,1], modes[inds])
coes.(non_measurable_modes)

coefficients = vcat(coefficients, [abs_coes.(non_measurable_modes) .* exp.(im .* non_measurable_modes) non_measurable_modes .* 0.0 .+ 0.0im])


# coefficients = vcat(coefficients, [coes.(non_measurable_modes) non_measurable_modes .* 0.0 .+ 0.0im])
modes = [modes; non_measurable_modes]

inds = sortperm_modes(modes)
modes = modes[inds] 
coefficients = coefficients[inds,:]

inds = sortperm(loading_profile.modes)
loading_profile.modes[inds]


fun = abs
fun = angle
fun = real

plot(loading_profile.modes[inds],loading_profile.coefficients[inds,1] .|> fun)
scatter!(modes[inds],coefficients[inds,1] .|> fun)
plot!(ylims = (-0.04,0.04), xlims = (10,70))

scatter(modes[inds],coefficients[inds,1] .|> fun)

scatter!(loading_profile.modes[inds],loading_profile2.coefficients[inds,1] .|> fun)
plot!(ylims = (-0.04,0.04), xlims = (10,70))

loading_profile3 = BoundaryData(bc1_forward, 
    θs = θs, 
    coefficients = coefficients,
    modes = modes
)

loading_profile3 = fouriermodes_to_fields(loading_profile3)


loading_profile2 = fouriermodes_to_fields(loading_profile2)

plot(loading_profile.θs, loading_profile.fields .|> abs)

plot!(loading_profile.θs, loading_profile2.fields .|> abs, linestyle = :dash)
plot!(loading_profile3.θs, loading_profile3.fields .|> abs, linestyle = :dash)

plot(loading_profile3.θs, loading_profile3.fields .|> imag, linestyle = :dash)

