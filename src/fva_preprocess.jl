# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)
function fva_preprocess(S, b, lb, ub, rxns; 
        check_obj = nothing,
        verbose = true, eps = 0.0, 
        batchlen = 50,
        check_obj_atol::Real = 1e-4,
        ignore = [], # rxns skip preprocess
        protect = [], # rxns skip blocking
    )

    # FVA
    _bidx = trues(size(S, 2))
    _bidx[ignore] .= false
    non_ignored = findall(_bidx)

    fvalb, fvaub = (lb, ub) .|> copy
    fvalb[non_ignored], fvaub[non_ignored] = fva(
        S, b, lb, ub, non_ignored; 
        check_obj_atol, check_obj, 
        batchlen, verbose
    )

    return del_blocked(S, b, fvalb, fvaub, rxns; eps, protect)
end


function fva_preprocess(metnet::MetNet; 
            check_obj = nothing,
            verbose = true, eps = 0.0, 
            return_blocked = false,
            batchlen = 50,
            check_obj_atol::Real = 1e-4,
            ignore = [], # rxns skip preprocess
            protect = [] # rxns skip blocking
        )
    ignore = map((r) -> rxnindex(metnet, r), ignore)
    protect = map((r) -> rxnindex(metnet, r), protect)
    check_obj = isnothing(check_obj) ? nothing : rxnindex(metnet, check_obj)
    
    model_fields = _extract_dense(metnet, [:S, :b, :lb, :ub])
    rxns = eachindex(metnet.rxns)
    S, b, lb, ub, rxnis, blocked = fva_preprocess(
        model_fields..., rxns; 
        check_obj, verbose, eps, 
        ignore, protect, batchlen, 
        check_obj_atol
    );
    
    rxns = metnet.rxns[rxnis]
    metnet = MetNet(metnet; S, b, lb, ub, rxns)
    return return_blocked ? (metnet, blocked) : metnet
end