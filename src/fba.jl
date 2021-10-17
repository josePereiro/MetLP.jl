function fba(S, b, lb, ub, obj_idx::Integer; 
        sense = MAX_SENSE, 
        sv = zeros(size(S, 2)), 
        drop_LPsol = true
    )
    sv[obj_idx] = sense
    LPsol = linprog(
        sv, # Opt sense vector 
        S, # Stoichiometric matrix
        b, # row lb
        b, # row ub
        lb, # column lb
        ub, # column ub
        ClpSolver()
    )
    v = LPsol.sol
    LPsol = drop_LPsol ? nothing : LPsol
    return isempty(v) ? 
        FBAOut(fill(NaN, size(S, 2)), NaN, obj_idx, LPsol) : 
        FBAOut(v, v[obj_idx], obj_idx, LPsol)
end

function fba(S, b, lb, ub, idx1::Integer, idx2::Integer;
        sense1 = MAX_SENSE,
        sense2 = MIN_SENSE,
        btol = 0.0, # tol of the fixation
        kwargs...
    )
    # maximizing obj
    FBAOut1 = fba(S, b, lb, ub, idx1; sense = sense1, kwargs...)
    obj_val = FBAOut1.obj_val
    # fix obj
    # TODO: do not copy here
    lb_ = copy(lb)
    ub_ = copy(ub)
    dflux = abs(obj_val * btol)
    lb_[idx1] = obj_val - dflux
    ub_[idx1] = obj_val + dflux
    # minimize cost
    return fba(S, b, lb_, ub_, idx2; sense = sense2, kwargs...)
end

function fba!(model::MetNet, obj_ider::IDER_TYPE; 
        sv = model.c, kwargs...
    )
    obj_idx = rxnindex(model, obj_ider)
    model_fields = _extract_dense(model, [:S, :b, :lb, :ub])
    return fba(model_fields..., obj_idx; sv, kwargs...)
end

function fba!(model::MetNet, ider1::IDER_TYPE, ider2::IDER_TYPE; 
        sv = model.c, 
        kwargs...
    )
    idx1 = rxnindex(model, ider1)
    idx2 = rxnindex(model, ider2)
    model_fields = _extract_dense(model, [:S, :b, :lb, :ub])
    return fba(model_fields..., idx1, idx2; sv, kwargs...)
end