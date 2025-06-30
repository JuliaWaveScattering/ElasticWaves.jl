
# Create plots of specific modes

using ElasticWaves
using Test, Statistics, LinearAlgebra, MultipleScattering
using Plots

medium = Elastic(2; ρ = 7000.0, cp = 5000.0 - 0.0im, cs = 3500.0 - 0.0im)

inner_radius = 1.0
bearing = RollerBearing(medium = medium, 
    inner_radius = inner_radius, outer_radius = inner_radius + 0.3, 
    number_of_rollers = 0,
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
    
    plot(result_pressure, result_pressure.ω[1];
        seriestype=:heatmap
        , field_apply = f -> real(f[1])
        # , leg = false,
    );
    plot!(
        # xlims = (-bearing.inner_radius * 0.4, bearing.inner_radius * 0.4), 
        # ylims = (bearing.inner_radius * 0.8, bearing.outer_radius * 1.1), 
        frame = :none, title="", xguide ="",yguide =""
    )

    dkr = 8.0;
    dn = 4
    
    dkrs = 2.5:dkr:(2dkr+2.5) |> collect
    dns = 2.0:dn:(dn+2) |> collect
    
    dns2 = reverse(Int[n for kr in dkrs, n in dns])
    dkrs2 = [kr for kr in dkrs, n in dns]


results = map(eachindex(dkrs2)) do i

    ω = real(medium.cp) * dkrs2[i];
    
    bd1 = BoundaryData(TractionBoundary(inner=true); coefficients =  [0.0 0.0], modes = [dns2[i]])
    bd2 = BoundaryData(TractionBoundary(outer=true); coefficients =  [1.0 1.0], modes = [dns2[i]])

    method = ModalMethod(
        tol = tol,
        only_stable_modes = false 
    )
    sim = BearingSimulation(ω, bearing, bd1, bd2; 
        method = method,
        nondimensionalise = true
    ) 
    wave = ElasticWave(sim)

    result_traction = field(wave, bearing, TractionType(); res = 250);
    fs = [f[1] for f in result_traction.field];
    result_pressure = FrequencySimulationResult(fs, result_traction.x, [ω]);
end

gr(size = (1200,800))

ps = map(results) do r
    plot(r, r.ω[1];
        seriestype=:heatmap
        , field_apply = f -> real(f[1])
        , leg = false,
    );
    plot!(
        # xlims = (-bearing.inner_radius * 0.4, bearing.inner_radius * 0.4), 
        # ylims = (bearing.inner_radius * 0.8, bearing.outer_radius * 1.1), 
        frame = :none, title="", xguide ="",yguide =""
    )
    plot!(bearing)
end;

gr(size = (1200,1200))
plot(ps..., layout = (2, 3))

# savefig("docs/images/modes-r1-1.0-r2-1p3.pdf")