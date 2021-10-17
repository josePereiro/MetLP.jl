
# given a model and a new set of bounds check
# how they will affect a given obj reaction
function check_newbounds(S, b, lb, ub, newlb, newub, 
        check_obj::Int, idxs = eachindex(lb); 
        check_obj_atol = 1e-4,
        verbose = true,
        batchlen = 50
    )

    M, N = size(S)
    T = eltype(S)
    nths = nthreads()
    ref_obj_val = fba(S, b, lb, ub, check_obj).obj_val

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
        new_obj_val = fba(S, b, wlb, wub, check_obj).obj_val
        if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)

            # reset batch
            wlb[batch] = @view lb[batch]
            wub[batch] = @view ub[batch]

            # If whole batch fail, test independent idxs
            for idx in batch
                
                # check both first 
                wlb[idx], wub[idx] = newlb[idx], newub[idx]
                new_obj_val = fba(S, b, wlb, wub, check_obj).obj_val

                if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)
                    wlb[idx], wub[idx] = lb[idx], ub[idx] # reset both
                    
                    # if obj_val changed I check each bound
                    for (wcol, newb, oldb) in [(wlb, newlb[idx], lb[idx]), 
                                                (wub, newub[idx], ub[idx])]
                        wcol[idx] = newb
                        new_obj_val = fba(S, b, wlb, wub, check_obj).obj_val
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


# given a model and a new set of bounds check
# # how they will affect a given obj reaction
# function check_newbounds(S, b, lb, ub, newlb, newub, 
#     check_obj::Int, idxs = eachindex(lb); 
#     check_obj_atol = 1e-4,
#     verbose = true,
#     batchlen = 50)

#     M, N = size(S)
#     T = eltype(S)
#     ref_obj_val = fba(S, b, lb, ub, check_obj).obj_val

#     # thread environment (avoid race)
#     wlb, wub = deepcopy.([lb, ub])

#     icount = length(idxs)
#     verbose && (prog = Progress(icount; desc = "Checking bounds  "))
        
#     # If whole batch fail, test independent idxs
#     for idx in idxs

#         # check both first (this use the heuristic that only a few rxns will affect the biomass)
#         wlb[idx], wub[idx] = newlb[idx], newub[idx]
#         new_obj_val = fba(S, b, wlb, wub, check_obj).obj_val

#         if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)
#             wlb[idx], wub[idx] = lb[idx], ub[idx] # reset both
            
#             # if obj_val changed I check each bound
#             for (wcol, newb, oldb) in [(wlb, newlb[idx], lb[idx]), 
#                                         (wub, newub[idx], ub[idx])]
#                 wcol[idx] = newb
#                 new_obj_val = fba(S, b, wlb, wub, check_obj).obj_val
#                 if !isapprox(new_obj_val, ref_obj_val; atol = check_obj_atol)
#                     wcol[idx] = oldb # reset
#                 end
#             end
#         end
#         verbose && next!(prog; showvalues = [
#                 (:ref_obj_val, ref_obj_val), 
#                 (:new_obj_val, new_obj_val)
#             ] 
#         )
#     end

#     return wlb, wub
# end