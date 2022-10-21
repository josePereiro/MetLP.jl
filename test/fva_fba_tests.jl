function test_fva_consistency()

    # model
    model = MetNets.ecoli_core_model();
    m, n = size(model)

    # parameters
    niter = 10
    robj = MetNets.ECOLI_MODEL_BIOMASS_IDER
    nths = clamp(nthreads() - 1, 1, 5)
    atol = 1e-5

    # Test fba
    fbaout0 = MetLP.fba(model, robj)
    objval0 = MetLP.av(model, fbaout0, robj)
    @show objval0
    @test objval0 > 0.0
    
    # Testing fva consistency checking and not checking an obj value
    for (testname, check_obj) in [
            ("No-Checking Obj, FVA  ", nothing), 
            ("Checking Obj, FVA  ", robj)
        ]
        
        # fva
        println("-"^60)
        @info("Testing FVA")
        prog = Progress(niter; desc = testname)
        ref_lbs, ref_ubs = MetLP.fva(model; check_obj, nths, verbose = false)
        for _ in 1:niter
            lbs, ubs = MetLP.fva(model; check_obj, nths, verbose = false)
            @test all(isapprox.(lbs, ref_lbs; atol))
            @test all(isapprox.(ubs, ref_ubs; atol))
            next!(prog)
        end
        finish!(prog)
        
        @info("Testing FBA")
        for (rxni, rxn) in enumerate(model.rxns)

            println("-"^60)
            @show rxn
            @show model.lb[rxni]
            @show model.ub[rxni]

            # max
            fbaout = MetLP.fba(model, rxn, 1; sense1 = MetLP.MIN_SENSE)
            @test MetLP.av(model, fbaout, rxn) == MetLP.av(fbaout, rxni) # test interface
            val = MetLP.av(model, fbaout, rxn)
            @show val
            @show ref_lbs[rxni]
            @test all(isapprox(ref_lbs[rxni], val; atol))

            # min
            fbaout = MetLP.fba(model, rxn, 1; sense1 = MetLP.MAX_SENSE)
            val = MetLP.av(model, fbaout, rxn)
            @test val == MetLP.av(fbaout, rxni)
            @show val
            @show ref_ubs[rxni]
            @test all(isapprox(ref_ubs[rxni], val; atol))

        end
    end
end
test_fva_consistency()