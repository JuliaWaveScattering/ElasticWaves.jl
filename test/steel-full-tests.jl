# Here we test 

# @testset "Steel bearing" begin

# Well almost steel!
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

min_condition_number = 1e-4
ωs = ωms
Cs = map(ωs) do ω
    conds = map(ns) do n
        A = boundarycondition_system(ω, bearing, bc1_forward, bc2_forward, n) 
        SM = diagm([(4.0) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
        A = A * SM
        A |> det |> abs
    end
    # plot(ns,conds, ylims = (0,0.001))

    n = findfirst(conds .< min_condition_number) - 1
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

plot(θs,qs) 
