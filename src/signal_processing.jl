function fields_to_fouriermodes(boundarydata::BoundaryData, basis_order::Int = round(floor(length(boundarydata.θs)/2 - 1/2)) |> Int)
    modes = fields_to_fouriermodes(boundarydata.θs, boundarydata.fields, basis_order)
    
    return BoundaryData(boundarydata.boundarytype;
        fourier_modes = modes,
        fields = boundarydata.fields,
        θs = boundarydata.θs
    )
end

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
