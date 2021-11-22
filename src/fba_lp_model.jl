## ------------------------------------------------------------------------------
const LPModel = JuMP.Model

## ------------------------------------------------------------------------------
const _RXNS_VAR_KEY = :rxnv
const _STOI_CON_KEY = :stoi_conv
const _LB_CON_KEY = :lb_conv
const _UB_CON_KEY = :ub_conv

_length(lp_model::LPModel) = length(_get_vars(lp_model))

function _up_con_rhs!(lp_model::LPModel, conv_key::Symbol, vals::Vector, cidxs = 1:_length(lp_model))
    @assert length(vals) == length(cidxs)
    conv = lp_model[conv_key]
    for (ci, val) in zip(cidxs, vals)
        JuMP.set_normalized_rhs(conv[ci], val)
    end
    return lp_model
end
_up_con_rhs!(lp_model::LPModel, conv_key::Symbol, val::Number, cidx::Integer) = _up_con_rhs!(lp_model, conv_key, [val], [cidx])

_get_con_rhs(lp_model::LPModel, conv_key::Symbol) = JuMP.normalized_rhs.(lp_model[conv_key])

function _delete!(lp_model::LPModel, con)
    try
        JuMP.delete(lp_model, con)
        catch err; !(err isa JuMP.MOI.InvalidIndex) && rethrow(err)
    end
    return lp_model
end
_delete!(lp_model::LPModel, sym::Symbol) = (_delete!(lp_model, lp_model[sym]); JuMP.unregister(lp_model, sym))

_get_vars(lp_model::LPModel) = lp_model[_RXNS_VAR_KEY]
_set_vars!(lp_model::LPModel, N) = (lp_model[_RXNS_VAR_KEY] = @JuMP.variable(lp_model, [1:N]); lp_model)

## ------------------------------------------------------------------------------
function build_lp_model(S, b, lb, ub, solver)
    lp_model = LPModel(solver)
    JuMP.set_silent(lp_model)
    _set_vars!(lp_model, size(S, 2))
    lb!(lp_model, lb)
    ub!(lp_model, ub)
    set_stoi_con!(lp_model, S, b)
    return lp_model
end
build_lp_model(net::MetNets.MetNet, solver) = build_lp_model(net.S, net.b, net.lb, net.ub, solver)
build_lp_model(net::MetNets.MetNet; solver = Clp.Optimizer) = build_lp_model(net, solver)

## ------------------------------------------------------------------------------
function set_start_value(lp_model::LPModel, vals)
    xs = _get_vars(lp_model)
    for (x, val) in zip(xs, vals)
        JuMP.set_start_value(x, val)
    end
    return lp_model
end

## ------------------------------------------------------------------------------
function set_stoi_con!(lp_model::LPModel, S, b)
    haskey(lp_model, _STOI_CON_KEY) && _delete!(lp_model, _STOI_CON_KEY)
    x = _get_vars(lp_model)
    lp_model[_STOI_CON_KEY] = @JuMP.constraint(lp_model, S * x .- b .== 0.0)
    return lp_model
end

## ------------------------------------------------------------------------------
# lb
import MetNets.lb
import MetNets.lb!
lb(lp_model::LPModel) = JuMP.normalized_rhs.(lp_model[_LB_CON_KEY])
lb(lp_model::LPModel, ridx::Int) = JuMP.normalized_rhs(lp_model[_LB_CON_KEY][ridx])
lb(lp_model::LPModel, ridxs) = JuMP.normalized_rhs.(lp_model[_LB_CON_KEY][ridxs])
function lb!(lp_model::LPModel, cidxs, lb)
    if haskey(lp_model, _LB_CON_KEY)
        _up_con_rhs!(lp_model, _LB_CON_KEY, lb, cidxs)
    else
        x = _get_vars(lp_model)
        lp_model[_LB_CON_KEY] = @JuMP.constraint(lp_model, x .>= lb, base_name = string(_LB_CON_KEY))
    end
    return lp_model
end
lb!(lp_model::LPModel, lb) = lb!(lp_model, 1:_length(lp_model), lb)

## ------------------------------------------------------------------------------
# ub
import MetNets.ub
import MetNets.ub!
ub(lp_model::LPModel) = JuMP.normalized_rhs.(lp_model[_UB_CON_KEY])
ub(lp_model::LPModel, ridx::Int) = JuMP.normalized_rhs(lp_model[_UB_CON_KEY][ridx])
ub(lp_model::LPModel, ridxs) = JuMP.normalized_rhs.(lp_model[_UB_CON_KEY][ridxs])
function ub!(lp_model::LPModel, cidxs, ub)
    if haskey(lp_model, _UB_CON_KEY)
        _up_con_rhs!(lp_model, _UB_CON_KEY, ub, cidxs)
    else
        x = _get_vars(lp_model)
        lp_model[_UB_CON_KEY] = @JuMP.constraint(lp_model, x .<= ub, base_name = string(_UB_CON_KEY))
    end
    return lp_model
end
ub!(lp_model::LPModel, ub) = ub!(lp_model, 1:_length(lp_model), ub)

## ------------------------------------------------------------------------------
function optimize!(lp_model::LPModel, obj_idx; sense = JuMP.MOI.MAX_SENSE)
    x = _get_vars(lp_model)
    @JuMP.objective(lp_model, sense, x[obj_idx])
    JuMP.optimize!(lp_model)
    return lp_model
end

## ------------------------------------------------------------------------------
function FBAOut(lp_model::LPModel, obj_idx; drop_LPsol = true)
    N = length(_get_vars(lp_model))
    status = JuMP.termination_status(lp_model)
    LPsol = drop_LPsol ? nothing : lp_model
    if status != JuMP.MOI.OPTIMAL
        return FBAOut(fill(NaN, N), NaN, obj_idx, LPsol)
    end
    x = _get_vars(lp_model)
    v = JuMP.value.(x)
    return FBAOut(v, v[obj_idx], obj_idx, LPsol)
end