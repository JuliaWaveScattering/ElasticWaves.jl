
# Here we test 

using ElasticWaves
using Test, Statistics, LinearAlgebra, MultipleScattering
using Plots, DataInterpolations

# used attenuation from Ultrasonic attenuation in pearlitic steel (2014)
medium = Elastic(2; ρ = 7000.0, cp = 5000.0 - 7.0im, cs = 3500.0 - 5.5im)
medium = Elastic(2; ρ = 1000.0, cp = 1000.0 - 1.2im, cs = 700.0 - 1.1im)

# need a higher rotation speed Ω for the forward and inverse problem to be well posed. In practice, we need only the inverse problem to be well posed.
Ω = 2pi * 700 / 60 # large wind turbines rotate at about 15 rpm
Z = 10 

inner_radius = 1.0
bearing = RollerBearing(medium = medium, 
    inner_radius = inner_radius, outer_radius = 1.1, 
    angular_speed = Ω,  
    rollers_inside = true,
    number_of_rollers = Z,
    roller_radius = inner_radius / (typeof(inner_radius)(1.0) + 2.5*Z / (2pi))
)

plot(bearing)

## Types of boundary conditions for the forward and inverse problem
    const bc1_forward = TractionBoundary(inner=true)
    const bc2_forward = TractionBoundary(outer=true)

    const bc1_inverse = DisplacementBoundary(outer=true)
    const bc2_inverse = TractionBoundary(outer=true)

    const max_condition_number = 1e7
    const tol = max_condition_number * eps(Float64)

## Let's check the resonant modes of this raceway
    
    frequency_order = 20
    ωms = natural_frequencies(bearing, frequency_order)
    ωs = LinRange(ωms[1], ωms[end], 2000)

    modes = 0:12
    resonances_vec = map(ωs) do ω
        conds = map(modes) do n
            nondim_bearing = nondimensionalise(bearing, ω)
            A = boundarycondition_system(ω, nondim_bearing, bc1_forward, bc2_forward, n) 
            SM = diagm([(4.0) / sum(abs.(A[:,j])) for j in 1:size(A,2)])
            A = A * SM
            cond(A)
            # 1 / abs(det(A))
        end
    end

    resonances_mat = hcat(resonances_vec...)
    heatmap(ωs, modes, resonances_mat, c = cgrad(:thermal, rev = true), 
        clims = (0,10000))

    map(ωms) do ω
        plot!([ω, ω], [0.0, modes[end]], c = :lightgreen, lab = "", linestyle = :dash)
    end       
    plot!()

    # tap test 

    θs = LinRange(-pi,pi, 50)[1:end-1]
    bd1_for = BoundaryData(bc1_forward, 
        θs = θs,
        fields = [ [1.0; zeros(Complex{Float64},  length(θs)-1) ]  zeros(Complex{Float64}, θs |> length)]
    )
    bd2_for = BoundaryData(bc2_forward, 
        θs = θs,
        fields = [zeros(Complex{Float64},  θs |> length) zeros(Complex{Float64}, θs |> length)]
    )

    modal_method = ModalMethod(tol = tol, only_stable_modes = true)
    ωs = LinRange(ωms[1], ωms[end], 1000)

    ω = ωs[4]
    accelerations =  map(ωs) do ω
        forward_sim = BearingSimulation(ω, modal_method, bearing, bd1_for, bd2_for);
        wave = ElasticWave(forward_sim);

        displacement_bd = BoundaryData(bc1_inverse, bearing.inner_radius, θs, wave)
        ω^2 * norm(displacement_bd.fields[:,1])
    end
    

    maximum(accelerations)
    minimum(accelerations)

    plot(ωs, accelerations, yaxis=:log)
    map(ωms) do ω
        plot!([ω, ω], [0.01, maximum(accelerations)], 
            c = :lightgreen, lab = "", 
            linestyle = :dash
            , yaxis=:log)
    end       
    plot!()

    plot(ωs[1:300], accelerations[1:300])
    plot(ωs, accelerations)
    
    frequency_order = 16
    ms = 1:frequency_order
    ωms = natural_frequencies(bearing, frequency_order) |> collect

    # plot!(ylims = (0,5))


## start with stribeck, then add defects by multiplying
    segment_number = 40
    θs = LinRange(-pi,pi, segment_number * bearing.number_of_rollers)[1:end-1]

    q0 = 6.0; ε = 0.5;
    q0s = q0 .* (1 .+ 0im .- 1/(2ε) .* (1 .- cos.(θs))) .^ (10.0 / 9.0)
    q0s = real.(q0s)
    q0s = [ ((q < 0) ? 0.0 : q) for q in q0s]

    θos = [-1.0];
    # θos = [-1.0, 0.1];
    is = [findmin(abs.(θo .- θs))[2] for θo in θos]

    qs_defect = θs .* 0.0 |> collect;
    qs_defect[is] = - q0s[is]

    qs = q0s + qs_defect 
    plot(θs,qs)

    loading_profile = BoundaryData(bc1_forward, 
        θs = θs, 
        fields = hcat(q0s .+ 0.0im,0.0 .* q0s)
    )
    loading_profile = fields_to_fouriermodes(loading_profile)

    loading_profile_defect = BoundaryData(bc1_forward, 
        θs = θs, 
        fields = hcat(qs_defect .+ 0.0im,0.0 .* qs_defect)
    )
    loading_profile_defect = fields_to_fouriermodes(loading_profile_defect)

## Below we calculate the displacement on the outer raceway
bd2_for = BoundaryData(bc2_forward, 
    θs = θs, 
    fields = [zeros(Complex{Float64},  θs |> length) zeros(Complex{Float64}, θs |> length)]
)

bd2_inverse = bd2_for # a little bit of an inverse crime

function outer_acceleration(ωs; 
        bearing = bearing,
        numberofsensors = 2,
        maximum_mode = 120,
        tol = tol,
        rel_error = 0.00,
        slip_amplitude = 0.0,
        defect_size = 0.2,
        loading_profile_defect = loading_profile_defect,
        loading_profile = loading_profile,
    )

    loading_profile_total = BoundaryData(bc1_forward, 
        θs = loading_profile.θs, 
        modes = loading_profile.modes, 
        fields = loading_profile.fields + defect_size .* loading_profile_defect.fields,
        coefficients = loading_profile.coefficients + defect_size .* loading_profile_defect.coefficients
    )

    ZΩ = bearing.number_of_rollers * bearing.angular_speed 

    modal_method = ModalMethod(tol = tol, only_stable_modes = true, maximum_mode = maximum_mode)
    θs_inv = LinRange(-pi,pi, numberofsensors + 1)[1:end-1]

    fields = map(ωs) do ω

        bd1_for = BoundaryData(ω, bearing, loading_profile_total);

        # add slip error
        θo = rand(-2.4:0.001:1.8)
        slip_amplitude = slip_amplitude * rand([-1.0, 1.0])
        p(θ) = θ + slip_amplitude * exp(-8 * (θ - θo)^2)

        # plot(θs,p.(θs))

        # force for slipped bearing
        f = LinearInterpolation(bd1_for.fields[:,1], bd1_for.θs; extrapolate = true)
        fs = abs.(f.(bd1_for.θs)) .* f.(p.(bd1_for.θs)) ./ abs.(f.(p.(bd1_for.θs)))

        # plot(θs, real.(fs - bd1_for.fields[:,1]))
        # plot(θs, real.(fs))
        # plot!(θs, real.(bd1_for.fields[:,1]), linestyle = :dash)

        # check if frequency is a harmonic of bearings
        fs = (round(ω / ZΩ) ≈ ω / ZΩ) ? fs : (fs - bd1_for.fields[:,1])

        amp = mean(abs.(bd1_for.fields))
        error = rel_error .* amp .* (rand(Complex{Float64},size(bd1_for.fields)...) .- 0.5 .- 0.5im)

        bd1_for = BoundaryData(bc1_forward, 
            θs = bd1_for.θs, 
            fields = [fs fs .* 0.0] + error
        )

        forward_sim = BearingSimulation(ω, modal_method, bearing, bd1_for, bd2_for);
        wave = ElasticWave(forward_sim);

        # bd1_for_recover = BoundaryData(bc1_forward, bearing.inner_radius, θs, wave)
        # plot(bd1_for_recover.θs, abs.(bd1_for_recover.fields[:,1]))
        
        bd1_inverse = BoundaryData(bc1_inverse, bearing.inner_radius, θs_inv, wave)
        bd1_inverse.fields[:] .* ω^2
    end

    return hcat(fields...) |> transpose |> collect
end

f0_mat = outer_acceleration(ωms; slip_amplitude = 0.00, defect_size = 0.0, rel_error = 0.0, numberofsensors = 2)
plot(ωms, 1e9 .* abs.(f0_mat[:,1]), linestyle = :dash)
scatter!(ωms, 0.0 .* ωms, lab = "")

dω = ωms[2] - ωms[1]
# ωs = dω:(dω/10):ωms[end]
ωs = ωms[5]:(dω/4):ωms[end]
# ωs = dω:(dω/2):ωms[end]
ωs |> length
rel_error = 0.0

f0_mat = outer_acceleration(ωs; slip_amplitude = 0.00, defect_size = 0.0, rel_error = rel_error, numberofsensors = 2)

plot(ωs, 1e9 .* abs.(f0_mat[:,1]), linestyle = :dash)
scatter!(ωms[5:end], 0.0 .* ωms, lab = "")

iterations = 20
f1_mats = [
    outer_acceleration(ωs; slip_amplitude = 0.01, defect_size = 0.0, rel_error = rel_error, numberofsensors = 2)
for i = 1:iterations]
f2_mats = [
    outer_acceleration(ωs; slip_amplitude = 0.0, defect_size = 0.1, rel_error = rel_error, numberofsensors = 2)
for i = 1:iterations]
f3_mats = [
    outer_acceleration(ωs; slip_amplitude = 0.01, defect_size = 0.1, rel_error = rel_error, numberofsensors = 2)
for i = 1:iterations]

plot(ωs, 1e9 .* abs.(f0_mat[:,1]), linestyle = :solid, lab = "baseline")
plot!(ωs, 1e9 .* abs.(mean(f2_mats)[:,1]), linestyle = :solid, lab = "defect")
plot!(ωs, 1e9 .* abs.(mean(f1_mats)[:,1]), lab = "slip", linestyle = :dash)
plot!(ωs, 1e9 .* abs.(mean(f3_mats)[:,1]), lab = "slip & defect", linestyle = :dashdot)
scatter!(ωms[5:end], 0.0 .* ωms, lab = "")

# plot(ωs, 1e9 .* abs.(f0_mat[:,1]), linestyle = :dash, lab = "baseline")
plot!(ωs, 1e9 .* abs.(f1_mats[1][:,1]), lab = "slip")
# plot!(ωs, 1e9 .* abs.(f2_mats[1][:,1]), lab = "slip & defect")
# scatter!(ωms[5:end], 0.0 .* ωms, lab = "")

decimate = 4;

ω2s = LinRange(ωs[1],ωs[1] + decimate * (ωs[end]-ωs[1]), decimate * (length(ωs) - 1) + 1)
ω2s[1:length(ωs)] - collect(ωs) |> norm

# ωs[1]:(ωs[2] - ωs[1]):(2ωs[end] + ωs[2] - 2ωs[1])
ts = ω_to_t(ω2s)

Fs_vec = [f0_mat[:,1], mean(f1_mats)[:,1], mean(f2_mats)[:,1], mean(f3_mats)[:,1]];
Fs_vec = [[f; zeros(length(ω2s) - length(f))] for f in Fs_vec];
fs_vec = [frequency_to_time(F, ω2s, ts) for F in Fs_vec]

plot(ts, fs_vec[1], xlims = (0,0.01))
plot!(ts, 0.12 .* cos.(10 .* ts .* Ω), linestyle = :dash)


using DSP
env_vec = [abs2.(hilbert(fs)) for fs in fs_vec]
plot!(env_vec[3])

Fenvs_vec = [time_to_frequency(env, ts, ω2s) for env in env_vec]


plot(ω2s, abs.(Fs_vec[1]))
plot!(ω2s, abs.(Fenvs_vec[1]), xlim = (ωs[1],ωs[end]))

plot(ω2s, abs.(Fs_vec[2]))
plot!(ω2s, abs.(Fenvs_vec[2]), xlim = (ωs[1],ωs[end]))

plot(ω2s, abs.(Fs_vec[3]))
plot!(ω2s, abs.(Fenvs_vec[3]), xlim = (ωs[1],ωs[end]))


plot(ω2s, abs.(Fs_vec[1]), xlim = (ωs[1],ωs[end]), linestyle=:dash, lab = "loading")
plot!(ω2s, abs.(Fs_vec[4]), xlim = (ωs[1],ωs[end]), linestyle=:dash, lab = "slip defect")
plot!(ω2s, abs.(Fenvs_vec[3]), xlim = (ωs[1],ωs[end]), lab = "env defect")
plot!(ω2s, abs.(Fenvs_vec[4]), xlim = (ωs[1],ωs[end]), ylims = (0,0.003), linestyle=:dash, lab = "env slip defect")


scatter!(ωms[5:end], 0.0 .* ωms, lab = "")
# ω0s = 0:(ωs[2]-ωs[1]):((length(ωs)-1) * (ωs[2]-ωs[1]))
# ts = ω_to_t(ω0s)
# f0_time = frequency_to_time(f0_mat[:,1], ω0s, ts)
# plot(f0_time, xlim = (1,6))


using DSP, Plots
x = 14:0.04:50
y0 = sin.(x)
y = cos.(2.3 .* x) .* y0 .+ cos.(2.3 .* x) .* sin.(1.2x)
y = cos.(2.3 .* x) .* y0 .* sin.(1.4x)

env = abs.(hilbert(y))
env2 = abs.(hilbert(env))

plot(x, abs.(y0))
plot!(x, y, linestyle=:dash)
plot!(x, env, color = :green)
plot!(x, env2, color = :red)
# plot!(x, -env, color = :green)

# using ProfileView
# @time outer_displacement(ωms; slip_amplitude = 0.0, defect_size = 0.0, rel_error = 0.0, numberofsensors = 2)

# @time outer_displacement(ωms)

# @profview outer_displacement(ωms; slip_amplitude = 0.0, defect_size = 0.0, rel_error = 0.0, numberofsensors = 2)

# f_mat2 = deepcopy(f_mat)




f_mat = outer_displacement(ωs; rel_error = 0.1, numberofsensors = 2, maximum_mode = 120)
plot!(ωs, abs.(f_mat[:,1]), linestyle = :dash)

# plot the fundamental frequencies and defect modes

# plot(ωs, abs.(f_mat[:,4]))

f = [f_mat[:,1]; f_mat[:,1] .* 0]

ωs = LinRange(dω, dω * length(f), length(f))
# ωs = LinRange(0, dω * length(f) - dω, length(f))

ts = ω_to_t(ωs)
F = frequency_to_time(f, ωs)

plot(ts, F)

modes = intersect(loading_predict.modes, loading_profile.modes)
loading_true = select_modes(loading_profile,modes)
loading_predict = select_modes(loading_predict,modes)

norm(loading_true.coefficients - loading_predict.coefficients) / norm(loading_true.coefficients)

loading_predict = fouriermodes_to_fields(loading_predict, θs)

plot(loading_predict.θs, loading_predict.fields[:,1] .|> real)
plot!(θs,qs_defect, linestyle = :dash)

norm(loading_predict.fields[:,1]- qs_defect) / norm(loading_predict.fields[:,1])

# Do a ribbon plot 

data = map(1:250) do i
    real.(predict_boundary(0.01,0:11).fields[:,1])
end

## redo the simulation 100 times each with a different error

# calculate the mean and standard deviation of y for each value of x
y_mat = hcat(data...) # simpler to make a big matrix
ys_mean = mean(y_mat, dims = 2)[:]
ys_std = std(y_mat, dims = 2)[:]

h = 200
gr(linewidth = 1.0, size = (1.6 * h, h ))

plot(θs, abs.(ys_mean); ribbon = ys_std, fillalpha=.4, linewidth = 2.0, lab = "predicted")
plot!(θs,qs, linestyle = :dash, lab = "exact")
plot!(xlab = "θ", ylab = "loading profile")

# savefig("docs/images/steel-bearing-defect-11-modes.pdf")

# Next we simulate the measured displacement on the outer raceway and test the signal processing methods typically used.





using DSP
θs = 0.0:0.04:(8pi)
f_fault = sin.(7 .* θs) .^ 3 
fs =  #(3 .+ 0.5 .* rand(θs |> length)) .* (cos.(17 .* θs)) + 
    (2 .+ 0.2 .* rand(θs |> length)) .* (cos.(15 .* θs)) + 
    f_fault
env = hilbert(fs)
plot(fs)
plot!(abs.(env))

moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

plot(fs)
plot!(moving_average(fs,15))

plot(abs.(env) .- mean(abs.(env)))
plot!(f_fault)

plot(moving_average(fs,15) .- mean(moving_average(fs,15)))
plot!(f_fault)