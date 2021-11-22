module MetLP

    using ProgressMeter
    using Base.Threads
    using SparseArrays

    import MetNets
    import MetNets: metindex, rxnindex, 
        MetNet, IDER_TYPE, AbstractMetState, 
        del_blocked, bounds, fixxing
    import JuMP, Clp

    include("FBAOut.jl")
    include("getters.jl")
    include("summary.jl")
    include("base.jl")
    include("const.jl")
    include("utils.jl")
    include("fba.jl")
    include("fba_lp_model.jl")
    include("fva.jl")
    include("fva_preprocess.jl")
    include("check_newbounds.jl")
    include("projection2D.jl")
    # include("yLP.jl")

end
