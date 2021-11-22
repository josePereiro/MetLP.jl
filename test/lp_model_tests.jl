function run_lp_model_tests()
    model = MetNets.toy_model()
    M, N = size(model)
    lp_model = MetLP.build_lp_model(model)

    # bounds are correct
    @test all(MetLP._get_con_rhs(lp_model, MetLP._LB_CON_KEY) .== model.lb)
    @test all(MetLP._get_con_rhs(lp_model, MetLP._UB_CON_KEY) .== model.ub)
    
    # test bound reset
    MetLP.set_lb_con!(lp_model, zeros(N))
    @test all(MetLP._get_con_rhs(lp_model, MetLP._LB_CON_KEY) .== 0.0)
    MetLP.set_ub_con!(lp_model, ones(N))
    @test all(MetLP._get_con_rhs(lp_model, MetLP._UB_CON_KEY) .== 1.0)

end
run_lp_model_tests();