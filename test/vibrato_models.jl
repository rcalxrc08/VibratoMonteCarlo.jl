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
        (name = "antithetic MC, MC integration", mc = MonteCarloConfiguration(Nsim, Nstep; sampling = FinancialMonteCarlo.AntitheticMC()), vb_mc = VibratoMonteCarloStandard(Nsim2)),
        (name = "Sobol MC, MC integration", mc = MonteCarloConfiguration(Nsim, Nstep; sampling = FinancialMonteCarlo.SobolMC()), vb_mc = VibratoMonteCarloStandard(Nsim2)),
    ]

    for (name, mc, vb_mc) in configs
        @testset "$name" begin
            is_analytic = vb_mc isa VibratoMonteCarloAnalytic
            # With MC inner integration (fresh inner samples per path, payoff baseline) values and deltas are within 2e-3,
            # gammas within 4.1e-2 (BS, serial).
            rtol_second = is_analytic ? 1e-2 : 1e-1
            rel_err_at(i, res, ref) = rel_err(value_and_derivatives(res)[i], value_and_derivatives(ref)[i])

            @testset "Kou, vibrato saltando" begin
                res = pricer(kou_model, rfCurve, VibratoSaltando(mc, vb_mc), eu_opt)
                test_value_and_derivatives(res, kou_ref; rtol = 1e-2, rtol_second = rtol_second)
                # Plain Float64 input must give the value of the AD run.
                res_v = pricer(kou(σ, S0_v, d), rfCurve, VibratoSaltando(mc, vb_mc), eu_opt)
                @test res_v ≈ value_and_derivatives(res)[1] rtol = 1e-12
            end

            @testset "Black Scholes, vibrato" begin
                res = pricer(bs_model, rfCurve, Vibrato(mc, vb_mc), eu_opt)
                test_value_and_derivatives(res, bs_ref; rtol = 1e-2, rtol_second = 1e-1)
            end

            @testset "Heston, vibrato" begin
                res = pricer(heston_model, rfCurve, Vibrato(mc, vb_mc), eu_opt)
                # The pathwise MC reference has no reliable second derivative.
                # Always serial: the antithetic Heston of FinancialMonteCarlo uses correlation -ρ.
                heston_ref = pricer(heston_model, rfCurve, MonteCarloConfiguration(Nsim, Nstep), eu_opt)
                @test rel_err_at(1, res, heston_ref) < 1e-2
                @test rel_err_at(2, res, heston_ref) < 1e-2
            end

            @testset "Variance Gamma, vibrato saltando" begin
                res = pricer(vg_model, rfCurve, VibratoSaltando(mc, vb_mc), eu_opt)
                test_value_and_derivatives(res, vg_ref; rtol = 1e-2, rtol_second = rtol_second)
            end
        end
    end
end

# PureLRM on Heston integrates the payoff against its marginal density at T, the FFT density of FinancialFFT (#51): the
# Gaussian PureLRM of the other Itô processes needs σ and a Lévy drift correction, which Heston has not. Relative errors
# against Lewis (value, delta, gamma) ≤ 4.4e-4 on calls with 200 nodes; the tolerance is about 3 times that. Binaries
# converge only O(dx) across the strike on the FFT trapezoid (#54), so they are not checked here.
@testset "PureLRM on Heston through the FFT density" begin
    rfCurve = ZeroRate(0.01)
    lewis = LewisMethod(400.0, 200000)
    S0 = hyper(100.0, 1.0, 1.0, 0.0)
    mc = MonteCarloConfiguration(10, 2)
    for K in (80.0, 100.0, 120.0)
        opt = EuropeanOption(1.0, K)
        test_value_and_derivatives(pricer(heston(0.2, S0, 0.01), rfCurve, PureLRM(mc, VibratoMonteCarloAnalytic()), opt), pricer(heston(0.2, S0, 0.01), rfCurve, lewis, opt); rtol = 1.5e-3)
    end
    # the derivatives over the parameters of the variance too: the gradient in (S0, vol of vol) against Lewis. The vol of
    # vol derivative lives in the tails of the density: 2.9e-3 off with the default ±7 standard deviations (2.7e-3 with
    # 800 nodes on them), 1.4e-6 on ±10.
    price(method, x) = pricer(heston(x[2], x[1], 0.01), rfCurve, method, EuropeanOption(1.0, 100.0))
    @test ForwardDiff.gradient(x -> price(PureLRM(mc, VibratoMonteCarloAnalytic(800, -10.0, 10.0)), x), [100.0, 0.2]) ≈ ForwardDiff.gradient(x -> price(lewis, x), [100.0, 0.2]) rtol = 1e-5
end
