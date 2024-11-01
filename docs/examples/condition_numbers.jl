# Create heatmaps of the condition number

using ElasticWaves
using Test, Statistics, LinearAlgebra, MultipleScattering
using Plots

medium = Elastic(2; ρ = 7000.0, cp = 5000.0 - 0.0im, cs = 3500.0 - 0.0im)

inner_radius = 1.0
bearing = RollerBearing(medium = medium, 
    inner_radius = inner_radius, outer_radius = inner_radius + 0.1, 
    angular_speed = Ω,  
    rollers_inside = true,
    number_of_rollers = Z,
    roller_radius = inner_radius / (typeof(inner_radius)(1.0) + 2.5*Z / (2pi))
)

krs = LinRange(0,50,400);
ωs = real(medium.cp) .* krs

ns = 0:50

# Types of boundary conditions for the forward and inverse problem
bc1_forward = TractionBoundary(inner=true)
bc2_forward = TractionBoundary(outer=true)

conds = [
    begin
        try 
            nondim_bearing = nondimensionalise(bearing, ω)
            A = boundarycondition_system(ω, nondim_bearing, bc1_forward, bc2_forward, n) 
            SM = diagm([(4.0) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
            A = A * SM
            cond(A)
        catch
            Inf
        end    
    end
for n in ns, ω in ωs]


using LaTeXStrings

# pyplot(size = (300,250))
gr(size = (300,250))
# Plots.PyPlot.rc("text", usetex="true")

heatmap(krs, ns, conds, clims = (0,25), 
    c = cgrad(:navia, rev = true)
    , xlab = L"k_p r_1"
    , ylab = L"n")

# savefig("docs/images/condition-forward-steel-r1-1.0-r2-1p5.pdf")

# plot!(10:40, (10:40) .* 0.0 .+ 20.)
# plot!(xlims = (30,31), ylims = (17,21), clims = (0,7525))


bc1_inverse = DisplacementBoundary(inner=true)
bc2_inverse = TractionBoundary(inner=true)

conds = [
    begin
        try 
            nondim_bearing = nondimensionalise(bearing, ω)
            A = boundarycondition_system(ω, nondim_bearing, bc1_inverse, bc2_inverse, n) 
            SM = diagm([(4.0) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
            A = A * SM
            cond(A)
        catch
            Inf
        end    
    end
for n in ns, ω in ωs]

heatmap(krs, ns, conds, clims = (0,25), 
    c = cgrad(:navia, rev = true)
    , xlab = L"k_p r_1"
    , ylab = L"n")

# savefig("docs/images/condition-inverse-steel.pdf")


bc1_inverse = DisplacementBoundary(inner=true)
bc2_inverse = DisplacementBoundary(outer=true)

conds = [
    begin
        try 
            nondim_bearing = nondimensionalise(bearing, ω)
            A = boundarycondition_system(ω, nondim_bearing, bc1_inverse, bc2_inverse, n) 
            SM = diagm([(4.0) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
            A = A * SM
            cond(A)
        catch
            Inf
        end    
    end
for n in ns, ω in ωs]

heatmap(krs, ns, conds, clims = (0,25), 
    c = cgrad(:navia, rev = true)
    , xlab = L"k_p r_1"
    , ylab = L"n")
