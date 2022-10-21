## ------------------------------------------------------------------------------
function optimize_MOMA_obj!(lp_model::LPModel, vr::AbstractVector; sense = MIN_SENSE)
    v = _get_vars(lp_model)
    JuMP.@objective(lp_model, sense, (v .- vr)' * (v .- vr))
    JuMP.optimize!(lp_model)
    return lp_model
end

## ------------------------------------------------------------------
function MOMA_fba(net::MetNet, vr::AbstractVector; 
        sense = MIN_SENSE,
        solver = Ipopt.Optimizer, 
        lp_model = nothing, 
        drop_LPsol = true
    )

    M, N = size(net)
    
    # prepare LP Model
    lp_model = isnothing(lp_model) ? build_lp_model(net; solver) : lp_model
    v = _get_vars(lp_model)

    # Optimize v^2
    optimize_MOMA_obj!(lp_model, vr; sense)
    
    # return
    LPsol = drop_LPsol ? nothing : lp_model
    !solution_found(lp_model) && return FBAOut(fill(NaN, N), NaN, -1.0, LPsol)
    return FBAOut(JuMP.value.(v), JuMP.objective_value(lp_model), -1.0, LPsol)

end