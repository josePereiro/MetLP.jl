function fva(S, b, lb, ub, ridxs = eachindex(lb); 
        check_obj::Union{Nothing, Int} = nothing, 
        check_obj_atol::Real = 1e-4,
        verbose = true, 
        batchlen::Int = 50,
        on_non_optimal_sol::Function = (idx, lp_model) -> error("FBA failed, non OPTIMAL solution returned!!!"),
        nths = nthreads(), 
        solver = Clp.Optimizer
    )

    ridxs = collect(ridxs)

    # bounds
    fvalb = lb[ridxs]
    fvaub = ub[ridxs]

    verbose && (prog = Progress(length(ridxs); desc = "Doing FVA (-t$nths)  "))
    ii_ch = Channel{Int}(nths) do ch_
        for ii in eachindex(ridxs)
            put!(ch_, ii)
            verbose && next!(prog)
        end
    end

    @threads for _ in 1:nths

        lp_model = build_lp_model(S, b, lb, ub, solver)
        x = _get_vars(lp_model)

        for ii in ii_ch
            ridx = ridxs[ii]

            # optimize max
            JuMP.@objective(lp_model, MIN_SENSE, x[ridx])
            JuMP.optimize!(lp_model)
            status = JuMP.termination_status(lp_model)
            if status == JuMP.MOI.OPTIMAL
                fvalb[ii] = JuMP.objective_value(lp_model)
                else; on_non_optimal_sol(ridx, lp_model)
            end

            # optimize min
            JuMP.@objective(lp_model, MAX_SENSE, x[ridx])
            JuMP.optimize!(lp_model)
            status = JuMP.termination_status(lp_model)
            if status == JuMP.MOI.OPTIMAL
                fvaub[ii] = max(JuMP.objective_value(lp_model), fvalb[ii])
                else; on_non_optimal_sol(ridx, lp_model)
            end
        end
    end
    verbose && finish!(prog)

    # Check bounds
    if !isnothing(check_obj)
        return MetLP.check_newbounds(S, b, lb, ub, fvalb, fvaub, 
            check_obj, ridxs; check_obj_atol, verbose, batchlen
        )
    end

    return fvalb, fvaub
end

function fva(model::MetNet, iders = eachindex(model.lb); 
        check_obj = nothing, kwargs...
    ) 
    obj_idx = isnothing(check_obj) ? nothing : rxnindex(model, check_obj)
    idxs = [rxnindex(model, idx) for idx in iders]
    model_fields = _extract_dense(model, [:S, :b, :lb, :ub])
    return fva(model_fields..., idxs; check_obj = obj_idx, kwargs...);
end

fva(model::MetNet, ider::IDER_TYPE; kwargs...) = fva(model, [ider]; kwargs...)