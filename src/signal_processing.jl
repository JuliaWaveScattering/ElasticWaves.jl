function fields_to_fouriermodes(boundarydata::AbstractBoundaryData, basis_order::Int = round(floor(length(boundarydata.θs)/2 - 1/2)) |> Int)
    
    modes = fields_to_fouriermodes(boundarydata.θs, boundarydata.fields, basis_order)
    
    # creates a copy of boundarydata
    @reset boundarydata.fourier_modes = modes
    
    return boundarydata
end

function fouriermodes_to_fields(boundarydata::AbstractBoundaryData)
    fields = fouriermodes_to_fields(boundarydata.θs, boundarydata.fourier_modes)
    
    @reset boundarydata.fields = fields

    return boundarydata
end


function fields_to_fouriermodes(θs::AbstractVector, fields::AbstractArray, basis_order::Int)

    if 2basis_order + 1 > length(θs)
        error("Can not calculate up to basis_order = $basis_order of the fourier modes from only $(length(θs)) field points. Either descrease basis_order or increase the number of points in fields")
    end

    exps = [
        exp(im * θ * m)
    for θ in θs, m = -basis_order:basis_order];

    fouriermodes = exps \ fields
end

function fouriermodes_to_fields(θs::AbstractVector, fouriermodes::AbstractArray)

    basis_order = basislength_to_basisorder(PhysicalMedium{2,1}, size(fouriermodes,1))

    exps = [
        exp(im * θ * m)
    for θ in θs, m = -basis_order:basis_order];

    fields = exps * fouriermodes
end

import LinearAlgebra: normalize!

function normalize!(bb::BoundaryBasis)
    for bd in bb.basis

        # calculate normalising factor
        n = if !isempty(bd.fields)

            # approximate the integral norm of the fields
            fs = (bd.fields[1:end-1,:] + circshift(bd.fields,-1)[1:end-1,:]) ./ 2
            dθs = circshift(bd.θs,-1)[1:end-1] - bd.θs[1:end-1]
            sum(abs2.(fs) .* dθs)
        elseif !isempty(bd.fourier_modes)
            2pi * norm(bd.fourier_modes)^2 
        else 
            return bb   
        end    

        bd.fourier_modes[:] = bd.fourier_modes ./ n
        bd.fields[:] = bd.fields ./ n
    end

    return bb
end