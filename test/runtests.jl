using MetLP
using MetLP.MetNets
using MetLP.ProgressMeter
using Test

@testset "MetLP.jl" begin
    include("fva_tests.jl")
    # include("yLP_tests.jl")
end
