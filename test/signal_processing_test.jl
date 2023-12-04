@testset "signal processing" begin

    # test invertability
    θs = 0.0:0.1:2pi
    basis_order = Int((length(θs) - 1)/2)

    f1 = rand(length(θs)) .+ im.* rand(length(θs))
    f2 = rand(length(θs)) .+ im.* rand(length(θs))
    fields = hcat(f1,f2)

    modes = fields_to_fouriermodes(θs,fields,basis_order)
    fields2 = fouriermodes_to_fields(θs,modes)

    @test maximum(abs.(fields - fields2)) < 1e-12

    # invertability is lost if basis_order is smaller
    basis_order = Int((length(θs) - 1)/2) - 1
    modes = fields_to_fouriermodes(θs,fields,basis_order)
    fields2 = fouriermodes_to_fields(θs,modes)

    @test maximum(abs.(fields - fields2)) < 0.3
    @test norm(fields - fields2) / norm(fields) < 0.2

    # add forcing with exact fourier modes
    # still invertable if basis_order is large enough to include the highest mode 17
    basis_order = 17

    f1 = rand() .* cos.(θs .* 3) + im .* rand() .* sin.(θs .* 14)
    f2 = rand() .* cos.(- θs .* 17) + im .* rand() .* cos.(θs .* 2)
    fields = hcat(f1,f2)

    modes = fields_to_fouriermodes(θs,fields,basis_order)
    fields2 = fouriermodes_to_fields(θs,modes)

    @test maximum(abs.(fields - fields2)) < 1e-12

    # but invertable lost if basis_order is lower
    basis_order = 16

    modes = fields_to_fouriermodes(θs,fields,basis_order)
    fields2 = fouriermodes_to_fields(θs,modes)

    @test maximum(abs.(fields - fields2)) / maximum(abs.(fields)) < 1.0
    @test norm(fields - fields2) / norm(fields) < 1.0

end
