# Create heatmaps of the condition number

using ElasticWaves
using Test, Statistics, LinearAlgebra, MultipleScattering
using Plots

medium = Elastic(2; ρ = 7000.0, cp = 5000.0 - 0.0im, cs = 3500.0 - 0.0im)

inner_radius = 1.0
bearing = RollerBearing(medium = medium, 
    inner_radius = inner_radius, 
    outer_radius = inner_radius + 0.3, 
)

krs = LinRange(0,50,400);
krs = LinRange(0,25,300);
# krs = LinRange(0,5,240);
ωs = real(medium.cp) .* krs

ns = 0:40
ns = 0:20

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
gr(size = (250,250))
gr(size = (310,250))
# Plots.PyPlot.rc("text", usetex="true")

heatmap(krs, ns, conds, clims = (0,25), 
    c = cgrad(:navia, rev = true)
    , xlab = L"k_p r_1"
    , ylab = L"n")

# savefig("docs/images/condition-forward-steel-r1-1.0-r2-1p5.pdf")

dkr = 8.0;
dn = 4

dkrs = 2.5:dkr:(2dkr+2.5) |> collect
dns = 2.0:dn:(dn+2) |> collect

dns2 = [n for kr in dkrs, n in dns]
dkrs2 = [kr for kr in dkrs, n in dns]

scatter!([kr for n in dns, kr in dkrs][:], [n for n in dns, kr in dkrs][:], lab = "")

# plot!(10:40, (10:40) .* 0.0 .+ 20.)
plot!(xlims = (0,20), ylims = (0,20), aspect_ratio = 1.0)
# savefig("docs/images/condition-scat-forward-steel-r1-1.0-r2-1p3.pdf")


bc1_inverse = DisplacementBoundary(outer=true)
bc2_inverse = TractionBoundary(outer=true)

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


## Check forward and inverse solve the same problem

krs = LinRange(0,50,500);
ωs = real(medium.cp) .* krs
ns = 0:50

n = -3;
ω = ωs[10]

δ = 1e-8
method = ModalMethod(only_stable_modes = false)

errors = [
    begin
        try 
            f0 = [1.0 1.0];
            amp = 0.02 .* mean(abs.(f0));
            error0 = amp .* ((rand(1,2) .- 0.5) + (rand(1,2) .- 0.5) .* im)
            error1 = amp .* ((rand(1,2) .- 0.5) + (rand(1,2) .- 0.5) .* im)
            error1 = error1 .* 0.0

            bd1 = BoundaryData(bc1_forward; coefficients =  f0 + error0, modes = [n])
            bd2 = BoundaryData(bc2_forward; coefficients =  [0.0 0.0] + error1, modes = [n])

            sim = BearingSimulation(ω, bearing, bd1, bd2; 
                method = method,
                nondimensionalise = true
            )
            wave = ElasticWave(sim);

            f1 = field_modes(wave, bearing.outer_radius, bc1_inverse.fieldtype)
            f2 = field_modes(wave, bearing.outer_radius, bc2_inverse.fieldtype)

            amp = 0.02 .* mean(abs.(f1));
            f1 = f1 + amp .* ((rand(1,2) .- 0.5) + (rand(1,2) .- 0.5) .* im)

            amp = 0.02 .* mean(abs.(f2));
            f2 = f2 + amp .* ((rand(1,2) .- 0.5) + (rand(1,2) .- 0.5) .* im)

            bd1 = BoundaryData(bc1_inverse; coefficients =  f1, modes = [n])
            bd2 = BoundaryData(bc2_inverse; coefficients =  f2, modes = [n])

            sim = BearingSimulation(ω, bearing, bd1, bd2; 
                method = method,
                nondimensionalise = true
            )
            wave = ElasticWave(sim);

            f1 = field_modes(wave, bearing.inner_radius, bc1_forward.fieldtype)
            f2 = field_modes(wave, bearing.outer_radius, bc2_forward.fieldtype)

            norm(f1 - f0) / sqrt(2)
            # norm(f2 - [0.0 0.0]) / sqrt(2)
        catch
            Inf
        end    
    end
for n in ns, ω in ωs]

# gr(size = (360,200))
gr(size = (310,250))

colorbar_ticks=(0:0.1:0.5, string.(round.(Int, (0:0.1:0.5) .* 100), "%"))

heatmap(krs, ns, errors, clims = (0,0.2), xlims = (0,50), ylims = (0,50)
    , colorbar_ticks = colorbar_ticks
    , aspect_ratio = 1.0
    , c = cgrad(:navia, rev = true)
    , xlab = L"k_p r_1"
    , ylab = L"n")

# savefig("docs/images/forward-inverse-cond-r1-1p0-r1-1p1.pdf")