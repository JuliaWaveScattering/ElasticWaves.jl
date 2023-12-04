@testset "traction and displacement formulas" begin

    medium = Elastic(2; ρ=7800.0, cp=5000.0, cs=3500.0)
    ω = 1.1
    basis_order = 4;

    basis_scaling = (1.0 ./  (abs.(-basis_order:basis_order) .+ 1).^2); 

    pcoes = (rand(2,2basis_order + 1) .- 0.5) + (rand(2,2basis_order + 1) .- 0.5) .* im
    pcoes = [pcoes[i] * basis_scaling[i[2]] for i in CartesianIndices(pcoes)]
    pressure = HelmholtzPotential{2}(medium.cp, ω / medium.cp, pcoes)
    
    scoes = (rand(2,2basis_order + 1) .- 0.5) + (rand(2,2basis_order + 1) .- 0.5) .* im
    scoes = [scoes[i] * basis_scaling[i[2]] for i in CartesianIndices(pcoes)]
   
    shear = HelmholtzPotential{2}(medium.cs, ω / medium.cs, scoes)

    wave = ElasticWave(ω, medium, [pressure, shear])

    # Let us implement an approximation for the displacement u = ∇φ + ∇x[0,0,ψ]
    # field(potential::HelmholtzPotential{Dim}, x::AbstractVector{T})

    
    function grad_pressure(pot::HelmholtzPotential, rθs::Vector; dr=1e-8, dθ=1e-8)

        fr1s = [field(pot, radial_to_cartesian_coordinates(rθs[i] - [dr, 0])) for i in eachindex(rθs)]
        fr2s = [field(pot, radial_to_cartesian_coordinates(rθs[i] + [dr, 0])) for i in eachindex(rθs)]
        
        fθ1s = [field(pot, radial_to_cartesian_coordinates(rθs[i] - [0, dθ])) for i in eachindex(rθs)]
        fθ2s = [field(pot, radial_to_cartesian_coordinates(rθs[i] + [0, dθ])) for i in eachindex(rθs)]

        dfdrs = (fr2s - fr1s) ./ (2dr)
        dfdθs = (fθ2s - fθ1s) ./ (2dθ)


        map(eachindex(rθs)) do i
            [dfdrs[i], dfdθs[i] / rθs[i][1]]
        end    
    end

    function curl_shear(pot::HelmholtzPotential, rθs::Vector; dr=1e-8, dθ=1e-8)

        fr1s = [field(pot, radial_to_cartesian_coordinates(rθs[i] - [dr, 0])) for i in eachindex(rθs)]
        fr2s = [field(pot, radial_to_cartesian_coordinates(rθs[i] + [dr, 0])) for i in eachindex(rθs)]

        fθ1s = [field(pot, radial_to_cartesian_coordinates(rθs[i] - [0, dθ])) for i in eachindex(rθs)]
        fθ2s = [field(pot, radial_to_cartesian_coordinates(rθs[i] + [0, dθ])) for i in eachindex(rθs)]

        dfdrs = (fr2s - fr1s) ./ (2dr)
        dfdθs = (fθ2s - fθ1s) ./ (2dθ)

        map(eachindex(rθs)) do i
            [dfdθs[i] / rθs[i][1], - dfdrs[i]]
        end    
    end

    function displacement_numerical(wave::ElasticWave, rθs; kws...)
        grad_pressure(wave.potentials[1], rθs; kws...) + curl_shear(wave.potentials[2], rθs; kws...)
    end    
    
    function strain_numerical(wave::ElasticWave, rθs; dr=1e-8, dθ=1e-8)
        us = displacement_numerical(wave, rθs, dr=dr, dθ=dθ)
        
        ur1s = displacement_numerical(wave, [rθ - [dr, 0] for rθ in rθs]; dr=dr, dθ=dθ)
        ur2s = displacement_numerical(wave, [rθ + [dr, 0] for rθ in rθs]; dr=dr, dθ=dθ)
        
        uθ1s = displacement_numerical(wave, [rθ - [0, dθ] for rθ in rθs]; dr=dr, dθ=dθ)
        uθ2s = displacement_numerical(wave, [rθ + [0, dθ] for rθ in rθs]; dr=dr, dθ=dθ)

        dudrs = (ur2s - ur1s) ./ (2dr)
        dudθs = (uθ2s - uθ1s) ./ (2dθ)
 
        return map(eachindex(rθs)) do i
            gradu = zeros(2,2) + zeros(2,2) .* im

            gradu[1,1] = dudrs[i][1]
            gradu[1,2] = (-us[i][2] + dudθs[i][1]) / rθs[i][1]
            gradu[2,1] = dudrs[i][2]
            gradu[2,2] = (us[i][1] + dudθs[i][2]) / rθs[i][1]

            gradu ./ 2 + transpose(gradu) ./2
        end    
    end    
    
    function stress_numerical(wave::ElasticWave, rθs; kws...)
        μ = wave.medium.ρ * wave.potentials[2].wavespeed^2;
        λ = wave.medium.ρ * wave.potentials[1].wavespeed^2 - 2μ;

        strains = strain_numerical(wave, rθs; kws...)
 
        return map(strains) do ε
            λ*tr(ε) .* I + 2μ .* ε
        end    
    end    

    rs = LinRange(600,700,100);
    θs = LinRange(0,2pi,100);

    rθs = [[rs[i], θs[i]] for i in eachindex(rs)]
    xs = radial_to_cartesian_coordinates.(rθs)

    dr = 1e-3
    dθ = 1e-5

    errors = norm.(displacement_numerical(wave, rθs; dr=dr, dθ=dθ) - [displacement(wave, x) for x in xs])
    max_error = maximum(errors) / mean(norm.([displacement(wave, x) for x in xs]));
    @test max_error < 1e-8


    σs = stress_numerical(wave, rθs; dr=dr, dθ=dθ);
    tractions = [σ[1:2, 1] for σ in σs];

    τs = [traction(wave, x) for x in xs];
    max_error = maximum(norm.(tractions - τs) / mean(norm.(tractions)))
    
    @test max_error < 1e-4
end  