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
            # With MC inner integration (fresh inner samples per path, payoff baseline) values and deltas are within 2e-3,
            # gammas within 4.1e-2 (BS, serial).
            rtol_second = is_analytic ? 1e-2 : 1e-1
            rel_err_at(i, res, ref) = rel_err(value_and_derivatives(res)[i], value_and_derivatives(ref)[i])

            @testset "Kou, vibrato saltando" begin
                res = vibrato_saltando(kou_model, rfCurve, mc, eu_opt, vb_mc)
                test_value_and_derivatives(res, kou_ref; rtol = 1e-2, rtol_second = rtol_second)
                # Plain Float64 input must give the value of the AD run.
                res_v = vibrato_saltando(kou(σ, S0_v, d), rfCurve, mc, eu_opt, vb_mc)
                @test res_v ≈ value_and_derivatives(res)[1] rtol = 1e-12
            end

            @testset "Black Scholes, vibrato" begin
                res = vibrato(bs_model, rfCurve, mc, eu_opt, vb_mc)
                test_value_and_derivatives(res, bs_ref; rtol = 1e-2, rtol_second = 1e-1)
            end

            @testset "Heston, vibrato" begin
                res = vibrato(heston_model, rfCurve, mc, eu_opt, vb_mc)
                # The pathwise MC reference has no reliable second derivative.
                # Always serial: the antithetic Heston of FinancialMonteCarlo uses correlation -ρ.
                heston_ref = pricer(heston_model, rfCurve, MonteCarloConfiguration(Nsim, Nstep), eu_opt)
                @test rel_err_at(1, res, heston_ref) < 1e-2
                @test rel_err_at(2, res, heston_ref) < 1e-2
            end

            @testset "Variance Gamma, vibrato saltando" begin
                res = vibrato_saltando(vg_model, rfCurve, mc, eu_opt, vb_mc)
                test_value_and_derivatives(res, vg_ref; rtol = 1e-2, rtol_second = rtol_second)
            end
        end
    end
end
