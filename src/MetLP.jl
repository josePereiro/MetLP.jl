module MetLP

    using ProgressMeter
    using Base.Threads
    using SparseArrays

    # TODO MathProgBase is deprecated, update the code!!!
    # Possibly make pull request for MetabolicEP too
    import MathProgBase
    import MathProgBase.HighLevelInterface: linprog
    import Clp: ClpSolver
    import MetNets
    import MetNets: metindex, rxnindex, 
        MetNet, IDER_TYPE, AbstractMetState, 
        del_blocked, bounds, fixxing
    # import JuMP, GLPK

    include("FBAOut.jl")
    include("base.jl")
    include("const.jl")
    include("utils.jl")
    include("fba.jl")
    include("fva.jl")
    include("fva_preprocess.jl")
    include("check_newbounds.jl")
    # include("yLP.jl")
    include("projection2D.jl")
    include("getters.jl")
    include("summary.jl")

end
