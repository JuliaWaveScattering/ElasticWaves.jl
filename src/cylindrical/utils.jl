function estimate_basisorder(ω::T, bearing::RollerBearing{T}; tol::T = 1e-5) where T

    @warn "this function is definitely wrong. I wrote it thinking about the inverse problem, but in general this needs to be re-thought and done more carefully. the basis_order should be increased until the solution has converged."

    kpa = bearing.outer_radius * ω / steel.cp
    ksa = bearing.outer_radius * ω / steel.cs
    ka = max(abs(kpa), abs(kpa))

    # estimate a maximum needed for most wave scattering
    max_basis_order = Int(round(5.0 * ka)) + 1

    # the inverse methods become ill posed when the bessel functions have changed significantly from one boundary to the other
    basis_orders = 0:max_basis_order
    ratios_hankel = [
        abs(hankelh1(m, ka * bearing.outer_radius) / hankelh1(m, ka * bearing.inner_radius))
    for m in basis_orders]

    ratios_bessel = [
        abs(besselj(m, ka * bearing.inner_radius) / besselj(m, ka * bearing.outer_radius))
    for m in basis_orders]

    i1 = findfirst(ratios_hankel .< tol)
    if isnothing(i1) i1 = max_basis_order end

    i2 = findfirst(ratios_bessel .< tol)
    if isnothing(i2) i2 = max_basis_order end

    basis_order = min(i1,i2)

    return basis_order
end
