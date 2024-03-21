@testset "signal processing" begin

    # test invertability
    θs = 0.0:0.1:2pi
    basis_order = Int((length(θs) - 1)/2)

    f1 = rand(length(θs)) .+ im.* rand(length(θs))
    f2 = rand(length(θs)) .+ im.* rand(length(θs))
    fields = hcat(f1,f2)

    modes = -basis_order:basis_order
    coefficients = fields_to_fouriermodes(θs,fields,modes)
    fields2 = fouriermodes_to_fields(θs,coefficients,modes)

    @test maximum(abs.(fields - fields2)) < 1e-12

    # invertability is lost if basis_order is smaller
    basis_order = Int((length(θs) - 1)/2) - 1
    modes = -basis_order:basis_order
    
    coefficients = fields_to_fouriermodes(θs,fields,modes)
    fields2 = fouriermodes_to_fields(θs,coefficients, modes)

    @test maximum(abs.(fields - fields2)) < 0.3
    @test norm(fields - fields2) / norm(fields) < 0.2

    # add forcing with exact fourier modes
    # still invertable if basis_order is large enough to include the highest mode 17
    basis_order = 17
    modes = -basis_order:basis_order

    f1 = rand() .* cos.(θs .* 3) + im .* rand() .* sin.(θs .* 14)
    f2 = rand() .* cos.(- θs .* 17) + im .* rand() .* cos.(θs .* 2)
    fields = hcat(f1,f2)

    coefficients = fields_to_fouriermodes(θs,fields,basis_order)
    fields2 = fouriermodes_to_fields(θs, coefficients, modes)

    @test maximum(abs.(fields - fields2)) < 1e-12

    # but invertable lost if basis_order is lower
    basis_order = 16

    modes = -basis_order:basis_order
    coefficients = fields_to_fouriermodes(θs,fields,modes)
    fields2 = fouriermodes_to_fields(θs,coefficients,modes)

    @test maximum(abs.(fields - fields2)) / maximum(abs.(fields)) < 1.2
    @test norm(fields - fields2) / norm(fields) < 1.0


    # Can invert between specific choice of fourier modes
    modes = [-2,4,6,7,8,14,15];
    coefficients = fields_to_fouriermodes(θs,fields,modes)
    
    # create a new fields with just the modes specificed
    fields = fouriermodes_to_fields(θs,coefficients,modes)
    
    # reproduce exactly the coefficients
    coefficients2 = fields_to_fouriermodes(θs,fields,modes)

    @test norm(coefficients2 - coefficients) / norm(coefficients) < 1e-12
    
    # also reproduce exactly the fields
    fields2 = fouriermodes_to_fields(θs,coefficients2,modes)
    
    @test norm(fields2 - fields) / norm(fields) < 1e-12
end
