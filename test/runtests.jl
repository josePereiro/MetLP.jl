using MetLP
using MetLP.MetNets
using MetLP.ProgressMeter
using Test
using Base.Threads

## ------------------------------------------------------
@testset "MetLP.jl" begin
    include("lp_model_tests.jl")
    include("fva_fba_tests.jl")
    # include("yLP_tests.jl")
end
