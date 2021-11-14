
# given a model and a new set of bounds check
# how they will affect a given obj reaction
function check_newbounds(S, b, lb, ub, newlb, newub, 
        check_obj::Int, idxs = eachindex(lb); 
        check_obj_atol = 1e-4,
        verbose = true,
        batchlen = 50, 
        solver = Clp.Optimizer
    )

    lp_model = build_lp_model(S, b, lb, ub, solver)
    fbaout = fba(lp_model, check_obj)
    isempty(fbaout) && error("fba returns is empty")
    ref_obj_val = objval(fbaout)

    # working copy
    wlb, wub = deepcopy.([lb, ub])

    icount = length(idxs)
    batchlen = max(1, min(batchlen, length(idxs)))
    batches = [idxs[i0:(min(i0 + batchlen - 1, icount))] for i0 in 1:batchlen:icount]
    verbose && (prog = Progress(length(batches); desc = "Checking bounds  "))
    for batch in batches
        
        # Test whole batch (this use the heuristic that only a few rxns will affect the biomass)
        wlb[batch] = @view newlb[batch]
        wub[batch] = @view newub[batch]

        # (TODO: This can be improved by updating only the relevant lp_model bounds)
        # Update bounds 
        set_lb_con!(lp_model, wlb)
        set_ub_con!(lp_model, wub)

        # fba
        fbaout = fba(lp_model, check_obj)
        new_obj_val = objval(fbaout)
        if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)

            # reset batch
            wlb[batch] = @view lb[batch]
            wub[batch] = @view ub[batch]

            # If whole batch fail, test independent idxs
            for idx in batch
                
                # check both first 
                wlb[idx], wub[idx] = newlb[idx], newub[idx]

                # Update bounds
                set_lb_con!(lp_model, wlb)
                set_ub_con!(lp_model, wub)

                # fba
                fbaout = fba(lp_model, check_obj)
                new_obj_val = objval(fbaout)

                if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)
                    wlb[idx], wub[idx] = lb[idx], ub[idx] # reset both
                    
                    # if obj_val changed I check each bound
                    for (wcol, newb, oldb) in [
                            (wlb, newlb[idx], lb[idx]), 
                            (wub, newub[idx], ub[idx])
                        ]
                        wcol[idx] = newb

                        # Update bounds
                        set_lb_con!(lp_model, wlb)
                        set_ub_con!(lp_model, wub)

                        # fba
                        fbaout = fba(lp_model, check_obj)
                        new_obj_val = objval(fbaout)
                        if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)
                            wcol[idx] = oldb # reset
                        end
                    end
                end
            end
        end

        verbose && next!(prog; showvalues = [
                (:ref_obj_val, ref_obj_val), 
                (:curr_obj_val, new_obj_val)
            ] 
        )

    end # for batch in batches

    verbose && finish!(prog)

    return wlb, wub
end
