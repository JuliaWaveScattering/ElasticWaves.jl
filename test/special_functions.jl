@testset "Special functions" begin

    ω = 5e4
    ω = 5e3

    steel = Elasticity(2; ρ = 7800.0, cp = 5000.0 - 0.1im, cs = 3500.0 - 0.1im)
    bearing = RollerBearing(medium=steel, r1=1.0, r2=2.0)

    # this non-dimensional number determines what basis_order is neeeded
    kpa = bearing.outer_radius * ω / steel.cp
    bearing.inner_radius * ω / steel.cp
    ksa = bearing.outer_radius * ω / steel.cs

    basis_order = Int(round(2.0 * max(abs(kpa),abs(ksa)))) + 1
    basis_length = basisorder_to_basislength(Acoustic{Float64,2}, basis_order)

    rs = LinRange(bearing.inner_radius,bearing.outer_radius,400);
    ms = 0:2:basis_order

    sqrts = 1 ./ sqrt.(rs .* ω / steel.cp)
    hs = [hankelh1.(m,rs .* ω / steel.cp) for m in ms]
    js = [besselj.(m,rs .* ω / steel.cp) for m in ms]
    hsx = [hankelh1x.(m,rs .* ω / steel.cp) for m in ms]

    using Plots

    plot(rs,abs.(hs[1]))
    plot!(rs,abs.(hs[2]))
    plot!(rs,abs.(hs[5]))
    plot!(rs,abs.(hs[15]))

    plot(rs,abs.(hsx[end]), xlims = (1.0,1.1))

    plot(rs,abs.(hs[4]), xlims = (1.0,1.5))
    plot(rs,abs.(js[end]))

    hankelh1.(basis_order,rs)


end
