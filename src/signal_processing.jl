function fields_to_fouriermodes(θs::AbstractVector, fields::AbstractArray, basis_order::Int)

    exps = [
        exp(im * θ * m)
    for θ in θs, m = -basis_order:basis_order];

    fouriermodes = exps \ fields
end

function fouriermodes_to_fields(θs::AbstractVector, fouriermodes::AbstractArray)

    basis_order = basislength_to_basisorder(Acoustic{Float64,2}, size(fouriermodes,1))

    exps = [
        exp(im * θ * m)
    for θ in θs, m = -basis_order:basis_order];

    fields = exps * fouriermodes
end
