function projection2D(model::MetNet, ider1::IDER_TYPE, ider2::IDER_TYPE; 
        l::Int = 20
    )
    
    ider1_L, ider1_U = bounds(model, ider1)
    ider1_range = range(ider1_L, ider1_U; length = l)
    ider2_L, ider2_U = bounds(model, ider2)
    ider2_range = range(ider2_L, ider2_U; length = l)

    proj = Dict{Float64, Vector{Float64}}()
    for flx1 in ider1_range
        try
            fixxing(model, ider1, flx1) do
                flx2_lb, flx2_ub = fva(model, ider2) .|> first
                proj[flx1] = [flx2_lb, flx2_ub]
            end
        catch err; end
    end
    return proj
end