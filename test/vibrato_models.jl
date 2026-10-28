@isdefined(TEST_UTILS_LOADED) || include("utils.jl")

@testset "vibrato on single name models" begin
    S0_v = 100.0
    S0 = hyper(S0_v, 1.0, 1.0, 0.0)
    K = 80.0
    r = 0.01
    T = 1.0
    d = 0.01
    σ = 0.2
    Nsim = 100_000
    Nstep = 300
    Nsim2 = 200
    rfCurve = ZeroRate(r)
    eu_opt = EuropeanOption(T, K)
    lewis = LewisMethod(400.0, 200000)

    kou_model = kou(σ, S0, d)
    bs_model = BlackScholesProcess(σ, Underlying(S0, d))
    heston_model = heston(σ, S0, d)
    vg_model = variance_gamma(σ, S0, d)
    kou_ref = pricer(kou_model, rfCurve, lewis, eu_opt)
    bs_ref = blsprice(S0, K, r, T, σ, d)
    vg_ref = pricer(vg_model, rfCurve, lewis, eu_opt)

    configs = [
        (name = "standard MC, analytic integration", mc = MonteCarloConfiguration(Nsim, Nstep), vb_mc = VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0)),
        (name = "standard MC, MC integration", mc = MonteCarloConfiguration(Nsim, Nstep), vb_mc = VibratoMonteCarloStandard(Nsim2)),
        (name = "antithetic MC, MC integration", mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC()), vb_mc = VibratoMonteCarloStandard(Nsim2)),
        (name = "Sobol MC, MC integration", mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SobolMode()), vb_mc = VibratoMonteCarloStandard(Nsim2)),
    ]

    for (name, mc, vb_mc) in configs
        @testset "$name" begin
            is_analytic = vb_mc isa VibratoMonteCarloAnalytic
            # With Sobol the inner integration points are all zero, so every derivative vanishes.
            is_sobol = mc.monteCarloMethod isa FinancialMonteCarlo.SobolMode
            delta_err(res, ref) = rel_err(value_and_derivatives(res)[2], value_and_derivatives(ref)[2])
            value_err(res, ref) = rel_err(value_and_derivatives(res)[1], value_and_derivatives(ref)[1])

            @testset "Kou, vibrato saltando" begin
                res = vibrato_saltando(kou_model, rfCurve, mc, eu_opt, vb_mc)
                if is_sobol
                    @test_broken value_err(res, kou_ref) < 1e-2
                else
                    @test value_err(res, kou_ref) < 1e-2
                end
                # Plain Float64 input must give the value of the AD run.
                res_v = vibrato_saltando(kou(σ, S0_v, d), rfCurve, mc, eu_opt, vb_mc)
                @test res_v ≈ value_and_derivatives(res)[1] rtol = 1e-12
                if is_analytic
                    test_value_and_derivatives(res, kou_ref; rtol = 1e-2)
                end
                is_sobol && @test_broken delta_err(res, kou_ref) < 5e-2
            end

            @testset "Black Scholes, vibrato" begin
                res = vibrato(bs_model, rfCurve, mc, eu_opt, vb_mc)
                @test value_err(res, bs_ref) < 1e-2
                if is_analytic
                    test_value_and_derivatives(res, bs_ref; rtol = 1e-2, rtol_second = 1e-1)
                end
                is_sobol && @test_broken delta_err(res, bs_ref) < 5e-2
            end

            @testset "Heston, vibrato" begin
                res = vibrato(heston_model, rfCurve, mc, eu_opt, vb_mc)
                # The pathwise MC reference has no reliable second derivative.
                heston_ref = pricer(heston_model, rfCurve, mc, eu_opt)
                @test value_err(res, heston_ref) < 1e-2
                is_analytic && @test delta_err(res, heston_ref) < 1e-2
                is_sobol && @test_broken delta_err(res, heston_ref) < 5e-2
            end

            @testset "Variance Gamma, vibrato saltando" begin
                res = vibrato_saltando(vg_model, rfCurve, mc, eu_opt, vb_mc)
                @test value_err(res, vg_ref) < 1e-2
                if is_analytic
                    @test delta_err(res, vg_ref) < 1e-2
                    @test_broken rel_err(value_and_derivatives(res)[3], value_and_derivatives(vg_ref)[3]) < 1e-1
                end
                is_sobol && @test_broken delta_err(res, vg_ref) < 5e-2
            end
        end
    end
end
