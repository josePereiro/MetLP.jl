function fba(S, b, lb, ub, obj_idx::Integer; 
        sense = MAX_SENSE, 
        on_empty_sol = () -> error("FBA failed, empty solution returned!!!")
    )
    sv = zeros(size(S, 2));
    sv[obj_idx] = sense
    sol = linprog(
        sv, # Opt sense vector 
        S, # Stoichiometric matrix
        b, # row lb
        b, # row ub
        lb, # column lb
        ub, # column ub
        ClpSolver());
    isempty(sol.sol) && on_empty_sol()
    return FBAOut(sol.sol, sol.sol[obj_idx], obj_idx, sol)
end

function fba(S, b, lb, ub, idx1::Integer, idx2::Integer;
        sense1 = MAX_SENSE, sense2 = MIN_SENSE,
        on_empty_sol = () -> error("FBA failed, empty solution returned!!!"), 
        btol = 0.0 # tol of the fixation
    )
    # maximizing obj
    FBAOut1 = fba(S, b, lb, ub, idx1; 
        sense = sense1, on_empty_sol
    )
    obj_val = FBAOut1.obj_val
    # fix obj
    lb_ = copy(lb)
    ub_ = copy(ub)
    dflux = abs(obj_val * btol)
    lb_[idx1] = obj_val - dflux
    ub_[idx1] = obj_val + dflux
    # minimize cost
    return fba(S, b, lb_, ub_, idx2; 
        sense = sense2, on_empty_sol
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

