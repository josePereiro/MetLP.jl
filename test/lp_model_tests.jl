function run_lp_model_tests()
    metnet = MetNets.toy_model()
    M, N = size(metnet)
    lp_model = MetLP.build_lp_model(metnet)

    # bounds are correct
    @show MetLP.lb(lp_model)
    @show MetLP.lb(metnet)
    @show MetLP.ub(lp_model)
    @show MetLP.ub(metnet)
    @test all(MetLP.lb(lp_model) .== MetLP.lb(metnet))
    @test all(MetLP.ub(lp_model) .== MetLP.ub(metnet))
    
    # test bound reset
    MetLP.lb!(lp_model, zeros(N))
    @test all(MetLP.lb(lp_model) .== 0.0)
    MetLP.ub!(lp_model, ones(N))
    @test all(MetLP.ub(lp_model) .== 1.0)

end
run_lp_model_tests();