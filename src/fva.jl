function fva(S, b, lb, ub, idxs = eachindex(lb); 
        check_obj::Union{Nothing, Int} = nothing, 
        check_obj_atol::Real = 1e-4,
        verbose = true, batchlen::Int = 50,
        zeroth::Real = 1e-10,
        on_empty_sol::Function = (idx, sense) -> error("FBA failed, empty solution returned!!!")
    )

    M, N = size(S)
    T = eltype(S)
    nths = nthreads()

    # Channels (avoid race)
    env_pool = Dict()
    for tid in 1:nths
        get!(env_pool, tid, 
            (
                tsv = zeros(T, N), 
                tfvalb = Dict{Int, T}(), 
                tfvaub = Dict{Int, T}()
            )
        )
    end

    icount = length(idxs)
    batchlen = max(1, min(batchlen, length(idxs)))
    batches = [idxs[i0:(min(i0 + batchlen - 1, icount))] for i0 in 1:batchlen:icount]
    verbose && (prog = Progress(length(batches); desc = "Doing FVA (-t$nths)  "))
    @threads for batch in batches
        
        # checks
        tid = threadid()
    
        # get environment
        tsv, tfvalb, tfvaub = env_pool[tid]
        
        verbose && next!(prog)
        for (fvacol, sense) in [(tfvalb, one(T)), (tfvaub, -one(T))]
            
            for idx in batch

                tsv[idx] = sense
                sol = linprog(
                    tsv, # Opt sense vector 
                    S, # Stoichiometric matrix
                    b, # row lb
                    b, # row ub
                    lb, # column lb
                    ub, # column ub
                    ClpSolver()
                )
                x = isempty(sol.sol) ? on_empty_sol(idx, sense) : sol.sol[idx]
                fvacol[idx] = abs(x) < zeroth ? zero(x) : x
                tsv[idx] = zero(sense)

            end
            
        end # for idx in idxs
    end
    verbose && finish!(prog)

    # Collect results
    merged_fvalb, merged_fvaub = Dict{Int, T}(), Dict{Int, T}()
    for (tid, (tsv, tfvalb, tfvaub)) in env_pool
        merge!(merged_fvalb, tfvalb)
        merge!(merged_fvaub, tfvaub)
    end

    # return just the indexed bounds
    fvalb, fvaub = getindex.([merged_fvalb], idxs), getindex.([merged_fvaub], idxs)

    # Check bounds
    if !isnothing(check_obj)
        return check_newbounds(S, b, lb, ub, fvalb, fvaub, 
            check_obj, idxs; check_obj_atol, verbose, batchlen
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