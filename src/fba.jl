function fba(lp_model::JuMP.Model, obj_idx::Integer; 
        sense = MAX_SENSE, drop_LPsol = true
    )
    
    optimize!(lp_model, obj_idx; sense)

    return FBAOut(lp_model, obj_idx; drop_LPsol)
end

function fba(S, b, lb, ub, obj_idx::Integer; 
        sense = MAX_SENSE, 
        drop_LPsol = true, 
        lp_model = nothing,
        up_lb_con = false,
        up_ub_con = false,
        up_stoi_con = false,
        solver = Clp.Optimizer
    )
    
    # build model
    if isnothing(lp_model)
        lp_model = build_lp_model(S, b, lb, ub, solver)
    else
        # update cons
        up_stoi_con && set_stoi_con!(lp_model, S, b)
        up_lb_con && lb!(lp_model, lb)
        up_ub_con && ub!(lp_model, ub)
    end

    fba(lp_model, obj_idx; sense, drop_LPsol)
end

function fba(lp_model::JuMP.Model, idx1::Integer, idx2::Integer; 
        sense1 = MAX_SENSE,
        sense2 = MIN_SENSE,
        btol = 0.0, # tol of the fixation
        drop_LPsol = true
    )

    # backup initial bounds
    lb0 = lb(lp_model, idx1)
    ub0 = ub(lp_model, idx1)

    # optimization idx1
    fbaout1 = fba(lp_model, idx1; sense = sense1, drop_LPsol)
    isempty(fbaout1) && return fbaout1
    val1 = objval(fbaout1)

    # fix obj1
    lb!(lp_model, idx1, val1 - btol)
    ub!(lp_model, idx1, val1 + btol)

    # optimization idx2
    fbaout2 = fba(lp_model, idx2; sense = sense2, drop_LPsol)

    # set back bounds
    lb!(lp_model, idx1, lb0)
    ub!(lp_model, idx1, ub0)

    return fbaout2
end

function fba(S, b, lb, ub, idx1::Integer, idx2::Integer;
        sense1 = MAX_SENSE,
        sense2 = MIN_SENSE,
        btol = 0.0, # tol of the fixation
        drop_LPsol = true,
        lp_model = nothing,
        up_lb_con = false,
        up_ub_con = false,
        up_stoi_con = false,
        solver = Clp.Optimizer
    )

    # build model
    if isnothing(lp_model)
        lp_model = build_lp_model(S, b, lb, ub, solver)
    else
        # update cons
        up_stoi_con && set_stoi_con!(lp_model, S, b)
        up_lb_con && lb!(lp_model, lb)
        up_ub_con && ub!(lp_model, ub)
    end
    
    return fba(lp_model, idx1, idx2; 
        sense1, sense2, btol, drop_LPsol
    )
    
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