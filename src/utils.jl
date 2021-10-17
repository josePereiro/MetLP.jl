function _extract_dense(model, fields)
    extracted = []
    for f in fields
        dat = getfield(model, f)
        if dat isa AbstractArray{<:Number} && issparse(dat)
            dat = Array(dat)
        end
        push!(extracted, dat)
    end
    return extracted
end