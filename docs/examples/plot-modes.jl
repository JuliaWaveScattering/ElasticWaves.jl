
# Create plots of specific modes

using ElasticWaves
using Test, Statistics, LinearAlgebra, MultipleScattering
using Plots

medium = Elastic(2; ρ = 7000.0, cp = 5000.0 - 0.0im, cs = 3500.0 - 0.0im)

inner_radius = 1.0
bearing = RollerBearing(medium = medium, 
    inner_radius = inner_radius, outer_radius = inner_radius * 1.1, 
    angular_speed = Ω,  
    rollers_inside = true,
    number_of_rollers = Z,
    roller_radius = inner_radius / (typeof(inner_radius)(1.0) + 2.5*Z / (2pi))
)

krs = LinRange(15,40,24400);
ωs = real(medium.cp) .* krs

# Types of boundary conditions for the forward and inverse problem
bc1_forward = TractionBoundary(inner=true)
bc2_forward = TractionBoundary(outer=true)

n = 20;
conds = [
    begin
        try 
            nondim_bearing = nondimensionalise(bearing, ω)
            A = boundarycondition_system(ω, nondim_bearing, bc1_forward, bc2_forward, n) 
            SM = diagm([(4.0) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
            A = A * SM
            cond(A)
        catch
            -1.0
        end    
    end
for ω in ωs]

i = findmax(conds)[2];
krs[i]
ω = ωs[i]

inds = findall(conds .> 0);
i = inds[findmin(conds[inds])[2]];

krs[i]
ω = ωs[i]

bd1 = BoundaryData(TractionBoundary(inner=true); coefficients =  [0.0 0.0], modes = [n])
bd2 = BoundaryData(TractionBoundary(outer=true); coefficients =  1e-4 .* [1.0 1.0], modes = [n])

    method = ModalMethod(only_stable_modes = false)
    sim = BearingSimulation(ω, bearing, bd1, bd2; 
        method = method
    ) 
    wave = ElasticWave(sim)

    result_traction = field(wave, bearing, TractionType(); res = 250);
    fs = [f[1] for f in result_traction.field];
    result_pressure = FrequencySimulationResult(fs, result_traction.x, [ω]);
    
    plot(result_pressure, ω;
        seriestype=:heatmap
        , field_apply = f -> real(f[1])
        # , leg = false,
    );
    plot!(
        xlims = (-bearing.inner_radius * 0.4, bearing.inner_radius * 0.4)
        , ylims = (bearing.inner_radius * 0.8, bearing.outer_radius * 1.1)
        ,frame = :none, title="", xguide ="",yguide =""
    )


sims = map(eachindex(ωs)) do i
    # the option only_stable_modes = false means the method will try to solve for modes which are ill posed 
    method = ModalMethod(regularisation_parameter = δs[i], 
        tol = tol, 
        only_stable_modes = false 
    )
    BearingSimulation(ωs[i], bearing, bd1, bd2; 
        method = method,
        nondimensionalise = true
    ) 
end

waves = [ElasticWave(s) for s in sims];

results = map(eachindex(ωs)) do i
    result_traction = field(inverse_waves[i], bearing, TractionType(); res = 250);
    fs = [f[1] for f in result_traction.field];
    result_pressure = FrequencySimulationResult(fs, result_traction.x, [ωs[i]]);
end