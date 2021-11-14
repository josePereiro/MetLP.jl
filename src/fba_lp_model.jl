const _RXNS_VAR_KEY = :rxnv
const _STOI_CON_KEY = :stoi_conv
const _LB_CON_KEY = :lb_conv
const _UB_CON_KEY = :ub_conv

_length(lp_model) = length(_get_vars(lp_model))

function _up_con_rhs!(lp_model::JuMP.Model, conv_key::Symbol, vals::Vector, cidxs = 1:_length(lp_model))
    @assert length(vals) == length(cidxs)
    conv = lp_model[conv_key]
    for (ci, val) in zip(cidxs, vals)
        JuMP.set_normalized_rhs(conv[ci], val)
    end
    return lp_model
end
_up_con_rhs!(lp_model::JuMP.Model, conv_key::Symbol, val::Number, cidx::Integer) = _up_con_rhs!(lp_model, conv_key, [val], [cidx])

function _delete!(lp_model, con)
    try
        JuMP.delete(lp_model, con)
        catch err; !(err isa JuMP.MOI.InvalidIndex) && rethrow(err)
    end
    return lp_model
end
_delete!(lp_model, sym::Symbol) = (_delete!(lp_model, lp_model[sym]); JuMP.unregister(lp_model, sym))

_get_vars(lp_model::JuMP.Model) = lp_model[_RXNS_VAR_KEY]
_set_vars!(lp_model::JuMP.Model, N) = (lp_model[_RXNS_VAR_KEY] = @JuMP.variable(lp_model, [1:N]); lp_model)

function build_lp_model(S, b, lb, ub, solver)
    lp_model = JuMP.Model(solver)
    JuMP.set_silent(lp_model)
    _set_vars!(lp_model, size(S, 2))
    set_lb_con!(lp_model, lb)
    set_ub_con!(lp_model, ub)
    set_stoi_con!(lp_model, S, b)
    return lp_model
end
build_lp_model(net::MetNets.MetNet, solver) = build_lp_model(net.S, net.b, net.lb, net.ub, solver) 

function set_stoi_con!(lp_model, S, b)
    haskey(lp_model, _STOI_CON_KEY) && _delete!(lp_model, _STOI_CON_KEY)
    x = _get_vars(lp_model)
    lp_model[_STOI_CON_KEY] = @JuMP.constraint(lp_model, S * x .- b .== 0.0)
    return lp_model
end

function set_lb_con!(lp_model, lb, cidxs = 1:_length(lp_model))
    if haskey(lp_model, _LB_CON_KEY)
        _up_con_rhs!(lp_model, _LB_CON_KEY, lb, cidxs)
    else
        x = _get_vars(lp_model)
        lp_model[_LB_CON_KEY] = @JuMP.constraint(lp_model, x .>= lb, base_name = string(_LB_CON_KEY))
    end
    return lp_model
end

function set_ub_con!(lp_model, ub, cidxs = 1:_length(lp_model))
    if haskey(lp_model, _UB_CON_KEY)
        _up_con_rhs!(lp_model, _UB_CON_KEY, ub, cidxs)
    else
        x = _get_vars(lp_model)
        lp_model[_UB_CON_KEY] = @JuMP.constraint(lp_model, x .<= ub, base_name = string(_UB_CON_KEY))
    end
    return lp_model
end

function _optimize!(lp_model, obj_idx; sense = JuMP.MOI.MAX_SENSE)
    x = _get_vars(lp_model)
    @JuMP.objective(lp_model, sense, x[obj_idx])
    JuMP.optimize!(lp_model)
    return lp_model
end

function FBAOut(lp_model::JuMP.Model, obj_idx; drop_LPsol = true)
    N = length(_get_vars(lp_model))
    LPsol = drop_LPsol ? nothing : lp_model
    status = JuMP.termination_status(lp_model)
    if status != JuMP.MOI.OPTIMAL
        return FBAOut(fill(NaN, N), NaN, obj_idx, LPsol)
    end
    x = _get_vars(lp_model)
    v = JuMP.value.(x)
    return FBAOut(v, v[obj_idx], obj_idx, LPsol)
end