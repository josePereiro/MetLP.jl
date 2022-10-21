## ------------------------------------------------------------------------------
function optimize_r2_fba_obj!(lp_model::LPModel; sense = MIN_SENSE)
    v = _get_vars(lp_model)
    JuMP.@objective(lp_model, sense, v' * v)
    JuMP.optimize!(lp_model)
    return lp_model
end

const _FIX_FBA_OBJ_CONST = :fix_fba

## ------------------------------------------------------------------
function r2_fba(net::MetNet, c::AbstractVector = net.c; 
        fba_sense = MAX_SENSE,
        r2_sense = MIN_SENSE,
        solver = Ipopt.Optimizer, 
        lp_model = nothing, 
        drop_LPsol = true
    )

    M, N = size(net)
    
    # prepare LP Model
    lp_model = isnothing(lp_model) ? build_lp_model(net; solver) : lp_model
    v = _get_vars(lp_model)
    haskey(lp_model, _FIX_FBA_OBJ_CONST) && _delete!(lp_model, _FIX_FBA_OBJ_CONST)

    # fba
    optimize_fba_obj!(lp_model, c; sense = fba_sense)
    LPsol = drop_LPsol ? nothing : lp_model
    !solution_found(lp_model) && return FBAOut(fill(NaN, N), NaN, findfirst(!iszero, c), LPsol)
    
    # fix objval
    objval = JuMP.objective_value(lp_model)
    lp_model[_FIX_FBA_OBJ_CONST] = JuMP.@constraint(lp_model, c' * v == objval)

    # Optimize v^2
    optimize_r2_fba_obj!(lp_model; sense = r2_sense)
    
    # return
    LPsol = drop_LPsol ? nothing : lp_model
    !solution_found(lp_model) && return FBAOut(fill(NaN, N), NaN, findfirst(!iszero, c), LPsol)
    return FBAOut(JuMP.value.(v), objval, findfirst(!iszero, c), LPsol)

end

## ------------------------------------------------------------------  
function r2_fba(net::MetNet, obj_ider::IDER_TYPE; kwargs...) 
    objidx = rxnindex(net, obj_ider)
    c = spzeros(size(net, 2))
    c[objidx] = 1.0
    return r2_fba(net, c; kwargs...)
end