function test_fva_consistency()

    # model
    model = MetNets.ecoli_core_model();
    m, n = size(model)

    # parameters
    niter = 100
    robj = MetNets.ECOLI_MODEL_BIOMASS_IDER
    nths = clamp(nthreads() - 1, 1, 5)
    atol = 1e-5
    
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

            # @show rxn

            # max
            fbaout = MetLP.fba(model, rxn; sense = MetLP.MIN_SENSE)
            val = MetLP.objval(fbaout)
            # @show val
            # @show ref_lbs[rxni]
            @test all(isapprox(ref_lbs[rxni], val; atol))

            # min
            fbaout = MetLP.fba(model, rxn; sense = MetLP.MAX_SENSE)
            val = MetLP.objval(fbaout)
            # @show val
            # @show ref_ubs[rxni]
            @test all(isapprox(ref_ubs[rxni], val; atol))

        end
    end

end
test_fva_consistency()