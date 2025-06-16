function fields_to_fouriermodes(boundarydata::AbstractBoundaryData, basis_order::Int = round(floor(length(boundarydata.θs)/2 - 1/2)) |> Int)

    fields_to_fouriermodes(boundarydata, -basis_order:basis_order)
end

function fields_to_fouriermodes(boundarydata::AbstractBoundaryData, modes::AbstractVector{Int})
    
    coefficients = fields_to_fouriermodes(boundarydata.θs, boundarydata.fields, modes)

    # creates a copy of boundarydata
    @reset boundarydata.coefficients = coefficients
    @reset boundarydata.modes = modes |> collect
    
    return boundarydata
end

function fouriermodes_to_fields(boundarydata::AbstractBoundaryData, θs::AbstractVector = boundarydata.θs)

    if θs |> isempty  
        θs = LinRange(0,2π, length(boundarydata.modes) + 1)[1:end-1]
    end    

    fields = fouriermodes_to_fields(θs, boundarydata.coefficients, boundarydata.modes)
    
    @reset boundarydata.fields = fields
    @reset boundarydata.θs = θs

    return boundarydata
end


function fields_to_fouriermodes(θs::AbstractVector, fields::AbstractArray, basis_order::Int)
    fields_to_fouriermodes(θs, fields, -basis_order:basis_order)
end

function fields_to_fouriermodes(θs::AbstractVector, fields::AbstractArray, modes::AbstractVector{Int})

    if length(modes) > length(θs)
        error("Can not calculate the modes = $modes of the Fourier series  from only $(length(θs)) field points. Either descrease the number of modes or increase the number of points in fields")
    end

    exps = [
        exp(im * θ * m)
    for θ in θs, m = modes];

    fouriermodes = exps \ fields
end

# the function below leads to unexpected errors
# function fouriermodes_to_fields(θs::AbstractVector, fouriermodes::AbstractArray, basis_order::Int = basislength_to_basisorder(PhysicalMedium{2,1}, size(fouriermodes,1)))

#     fouriermodes_to_fields(θs, fouriermodes, -basis_order:basis_order)
# end

function fouriermodes_to_fields(θs::AbstractVector, fouriermodes::AbstractArray, modes::AbstractVector{Int})

    exps = [
        exp(im * θ * m)
    for θ in θs, m = modes];

    fields = exps * fouriermodes
end

import LinearAlgebra: normalize!

"""
    sortperm_modes(modes::AbstractVector{Int})
    
Returns the indices `is` such that `modes[is]` is the standard ordering of the modes. Currently the lowest order modes are first.    
"""
function sortperm_modes(modes::AbstractVector{Int}) 

    is1 = sortperm(modes)
    is2 = sortperm(modes[is1], by = abs)

    return is1[is2]
end
is_standard_order(modes::AbstractVector{Int}) = (sortperm_modes(modes) == 1:length(modes))

function normalize!(bb::BoundaryBasis)
    for bd in bb.basis

        # calculate normalising factor
        n = if !isempty(bd.fields)

            # approximate the integral norm of the fields
            fs = (bd.fields[1:end-1,:] + circshift(bd.fields,-1)[1:end-1,:]) ./ 2
            dθs = circshift(bd.θs,-1)[1:end-1] - bd.θs[1:end-1]
            sum(abs2.(fs) .* dθs)
        elseif !isempty(bd.coefficients)
            2pi * norm(bd.coefficients)^2 
        else 
            return bb   
        end    

        bd.coefficients[:] = bd.coefficients ./ n
        bd.fields[:] = bd.fields ./ n
    end

    return bb
end