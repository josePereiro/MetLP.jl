function fba(S, b, lb, ub, obj_idx::Integer; 
        sense = MAX_SENSE, 
        drop_LPsol = true, 
        lp_model = nothing,
        up_lb_con = false,
        up_ub_con = false,
        up_stoi_con = false,
        solver = Clp.Optimizer
    )

    T = eltype(S)
    M, N = size(S)
    
    # build model
    if isnothing(lp_model)
        lp_model = _build_model(S, b, lb, ub, solver)
    end

    # update cons
    up_stoi_con && set_stoi_con!(lp_model, S, b)
    up_lb_con && set_lb_con!(lp_model, lb)
    up_ub_con && set_ub_con!(lp_model, ub)

    # setup
    JuMP.set_silent(lp_model)
    
    _optimize!(lp_model, obj_idx; sense)

    return FBAOut(lp_model, obj_idx; drop_LPsol)
end

function fba(S, b, lb, ub, idx1::Integer, idx2::Integer;
        sense1 = MAX_SENSE,
        sense2 = MIN_SENSE,
        btol = 0.0, # tol of the fixation
        drop_LPsol = true,
        lp_model = nothing,
        solver = Clp.Optimizer
    )

    # maximizing obj
    M, N = size(S)

    if isnothing(lp_model)
        lp_model = JuMP.Model(solver)
        @JuMP.variable(lp_model, lb[i] <= x[i=1:N] <= ub[i])
        @JuMP.constraint(lp_model, S * x .- b .== 0.0)
    end
    
    # setup
    x = lp_model[:x]
    LPsol = drop_LPsol ? nothing : lp_model
    JuMP.set_silent(lp_model)
    
    # optimization idx1
    @JuMP.objective(lp_model, sense1, x[idx1])
    JuMP.optimize!(lp_model)
    status = JuMP.termination_status(lp_model)
    if status != JuMP.MOI.OPTIMAL
        return FBAOut(fill(NaN, size(S, 2)), NaN, idx1, LPsol)
    end
    val1 = JuMP.objective_value(lp_model)
    
    # optimization idx2
    dflux = abs(val1 * btol)
    @JuMP.constraint(lp_model, val1 - dflux <= x[idx1] <= val1 + dflux)
    @JuMP.objective(lp_model, sense2, x[idx2])
    JuMP.optimize!(lp_model)
    status = JuMP.termination_status(lp_model)
    if status != JuMP.MOI.OPTIMAL
        return FBAOut(fill(NaN, size(S, 2)), NaN, idx1, LPsol)
    end
    
    v = JuMP.value.(x)
    FBAOut(v, v[idx1], idx1, LPsol)
end

function fba(model::MetNet, obj_ider::IDER_TYPE; kwargs...)
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