const _STOI_CON_KEY = :stoi_conv
const _LB_CON_KEY = :lb_conv
const _UB_CON_KEY = :ub_conv

_length(lp_model) = length(_get_var(lp_model))

function _up_con_rhs!(lp_model::JuMP.Model, conv_key::Symbol, vals::Vector, cidxs = 1:_length(lp_model))
    conv = lp_model[conv_key]
    for (vi, ci) in enumerate(cidxs)
        set_normalized_rhs(conv[ci], vals[vi])
    end
    return lp_model
end
_up_con_rhs!(lp_model::JuMP.Model, conv_key::Symbol, val, cidx::Integer) = _up_con_rhs!(lp_model, conv_key, [val], cidx)

function _delete!(lp_model, con)
    try
        JuMP.delete(lp_model, con)
        catch err; !(err isa JuMP.MOI.InvalidIndex) && rethrow(err)
    end
    return lp_model
end
_delete!(lp_model, sym::Symbol) = 
    haskey(lp_model, sym) ? 
        (_delete!(lp_model, lp_model[sym]); JuMP.unregister(model, sym)) : 
        lp_model

_get_var(lp_model::JuMP.Model) = lp_model[:x]

function _build_model(S, b, lb, ub, solver)
    M, N = size(S)
    lp_model = JuMP.Model(solver)
    JuMP.set_silent(lp_model)
    @JuMP.variable(lp_model, x[1:N])
    @JuMP.constraint(lp_model, lb_con, x .>= lb)
    @JuMP.constraint(lp_model, ub_con, x .<= ub)
    @JuMP.constraint(lp_model, stoi_con, S * x .- b .== 0.0)
    return lp_model
end

# TODO delete previous
function set_stoi_con!(lp_model, S, b)
    x = _get_var(lp_model)
    _delete!(lp_model, _STOI_CON_KEY)
    lp_model[_STOI_CON_KEY] = @JuMP.constraint(lp_model, S * x .- b .== 0.0)
    return lp_model
end

function set_lb_con!(lp_model, lb, cidxs = 1:_length(lp_model))
    x = _get_var(lp_model)
    if haskey(lp_model, _LB_CON_KEY)
        _up_con_rhs!(lp_model, _LB_CON_KEY, lb, cidxs)
    else
        lp_model[_LB_CON_KEY] = @JuMP.constraint(lp_model, x .<= lb, base_name = string(_LB_CON_KEY))
    end
    return lp_model
end

function set_ub_con!(lp_model, ub, cidxs = 1:_length(lp_model))
    x = _get_var(lp_model)
    if haskey(lp_model, _UB_CON_KEY)
        _up_con_rhs!(lp_model, _UB_CON_KEY, ub, cidxs)
    else
        lp_model[_UB_CON_KEY] = @JuMP.constraint(lp_model, x .<= ub, base_name = string(_UB_CON_KEY))
    end
    return lp_model
end

function _optimize!(lp_model, obj_idx; sense = JuMP.MOI.MAX_SENSE)
    x = _get_var(lp_model)
    @JuMP.objective(lp_model, sense, x[obj_idx])
    JuMP.optimize!(lp_model)
    return lp_model
end

function FBAOut(lp_model::JuMP.Model, obj_idx; drop_LPsol = true)
    N = length(_get_var(lp_model))
    LPsol = drop_LPsol ? nothing : lp_model
    status = JuMP.termination_status(lp_model)
    if status != JuMP.MOI.OPTIMAL
        return FBAOut(fill(NaN, N), NaN, obj_idx, LPsol)
    end
    v = JuMP.value.(lp_model[:x])
    return FBAOut(v, v[obj_idx], obj_idx, LPsol)
end