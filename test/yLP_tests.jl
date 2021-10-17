function test_yLP()
    
    model = Chemostat.Test.toy_model()
    glcidx = Chemostat.Utils.rxnindex(model, "gt")
    biomidx = Chemostat.Utils.rxnindex(model, "biom")
    
    d = zeros(length(model.c)); d[glcidx] = 1.0 
    model.c[biomidx] = 1.0
    
    status, yflxm, yield = Chemostat.LP.yLP(model, d);

    @test all(isapprox.(model.S * yflxm, model.b; atol = 1e-4))
    @test yield > 0.0
    @test status == Chemostat.LP.JuMP.MOI.OPTIMAL
end
test_yLP()