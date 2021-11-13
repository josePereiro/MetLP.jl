function fba(S, b, lb, ub, obj_idx::Integer; 
        sense = MAX_SENSE, 
        drop_LPsol = true, 
        solver = Clp.Optimizer, 
        on_non_optimal_sol::Function = (idx, lp_model) -> error("FBA failed, non OPTIMAL solution returned!!!"),
    )

    T = eltype(S)
    M, N = size(S)

    lp_model = JuMP.Model(solver)
    JuMP.set_silent(lp_model)
    @JuMP.variable(lp_model, lb[i] <= x[i=1:N] <= ub[i])
    @JuMP.constraint(lp_model, S * x .- b .== 0.0)
    @JuMP.objective(lp_model, sense, x[obj_idx])
    JuMP.optimize!(lp_model)
    status = JuMP.termination_status(lp_model)
    if status == JuMP.MOI.OPTIMAL
        v = JuMP.value.(x)
    else
        on_non_optimal_sol(obj_idx, lp_model)
        return FBAOut(T)
    end

    LPsol = drop_LPsol ? nothing : lp_model
    return isempty(v) ? 
        FBAOut(fill(NaN, size(S, 2)), NaN, obj_idx, LPsol) : 
        FBAOut(v, v[obj_idx], obj_idx, LPsol)
end

function fba(S, b, lb, ub, idx1::Integer, idx2::Integer;
        sense1 = MAX_SENSE,
        sense2 = MIN_SENSE,
        btol = 0.0, # tol of the fixation
        drop_LPsol = true,
        solver = Clp.Optimizer, 
        on_non_optimal_sol::Function = (idx, lp_model) -> error("FBA failed, non OPTIMAL solution returned!!!"),
    )

    # maximizing obj
    M, N = size(S)

    lp_model = JuMP.Model(solver)
    JuMP.set_silent(lp_model)
    @JuMP.variable(lp_model, lb[i] <= x[i=1:N] <= ub[i])
    @JuMP.constraint(lp_model, S * x .- b .== 0.0)
    
    # optimization idx1
    @JuMP.objective(lp_model, sense1, x[idx1])
    JuMP.optimize!(lp_model)
    status = JuMP.termination_status(lp_model)
    if status != JuMP.MOI.OPTIMAL
        on_non_optimal_sol(idx1, lp_model)
        return FBAOut(T)
    end
    val1 = JuMP.objective_value(lp_model)
    
    # optimization idx2
    dflux = abs(val1 * btol)
    @JuMP.constraint(lp_model, val1 - dflux <= x[idx1] <= val1 + dflux)
    @JuMP.objective(lp_model, sense2, x[idx2])
    JuMP.optimize!(lp_model)
    status = JuMP.termination_status(lp_model)
    if status != JuMP.MOI.OPTIMAL
        on_non_optimal_sol(idx2, lp_model)
        return FBAOut(T)
    end
    v = JuMP.value.(x)

    LPsol = drop_LPsol ? nothing : lp_model
    return isempty(v) ? 
        FBAOut(fill(NaN, size(S, 2)), NaN, idx1, LPsol) : 
        FBAOut(v, v[idx1], idx1, LPsol)
end

function fba(model::MetNet, obj_ider::IDER_TYPE; kwargs...
    )
    obj_idx = rxnindex(model, obj_ider)
    model_fields = _extract_dense(model, [:S, :b, :lb, :ub])
    return fba(model_fields..., obj_idx; kwargs...)
end

function fba(model::MetNet, ider1::IDER_TYPE, ider2::IDER_TYPE; kwargs...)
    idx1 = rxnindex(model, ider1)
    idx2 = rxnindex(model, ider2)
    model_fields = _extract_dense(model, [:S, :b, :lb, :ub])
    return fba(model_fields..., idx1, idx2; kwargs...)
end