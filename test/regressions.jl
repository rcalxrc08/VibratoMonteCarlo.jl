@isdefined(TEST_UTILS_LOADED) || include("utils.jl")
using Statistics

# Each case failed before the corresponding fix, with the error quoted in the comment.
@testset "regressions" begin
    r = 0.02
    eu_opt = EuropeanOption(1.0, 100.0)
    rfCurve = ZeroRate(r)

    # Antithetic paths used correlation -ρ between price and variance (corr ≈ +0.78 here).
    @testset "Heston antithetic correlation sign" begin
        h = HestonProcess(0.5, 0.2, 0.0, 2.0, -0.9, 0.04, Underlying(100.0, 0.0))
        mc = MonteCarloConfiguration(20_000, 50; sampling = FinancialMonteCarlo.AntitheticMC())
        v, S = VibratoMonteCarlo.simulate_with_last_vol(h, ZeroRate(0.0), mc, 1.0)
        @test cor(v, log.(S[:, end])) < -0.5
    end

    # The Analytic grid was in absolute log-price, the default one truncated at S ≈ 148 (27% too low here).
    @testset "PureLRM with default Analytic grid on Kou" begin
        model = kou(0.2, 100.0, 0.0)
        ref = pricer(model, rfCurve, LewisMethod(400.0, 200000), eu_opt)
        res = pricer(model, rfCurve, PureLRM(MonteCarloConfiguration(1_000, 10), VibratoMonteCarloAnalytic(2000)), eu_opt)
        @test rel_err(res, ref) < 1e-3
    end

    # The last jump time was not aligned with the grid point holding the jump, the drift was counted twice (+7.6% here).
    @testset "vibrato saltando on Kou with coarse grid" begin
        model = kou(0.2, 100.0, 0.0; λ = 50.0)
        ref = pricer(model, rfCurve, LewisMethod(600.0, 20000), eu_opt)
        mc = MonteCarloConfiguration(40_000, 10; sampling = FinancialMonteCarlo.AntitheticMC())
        res = pricer(model, rfCurve, VibratoSaltando(mc, VibratoMonteCarloAnalytic(200, -6.0, 6.0)), eu_opt)
        # Errors observed over 8 seeds ≤ 9.4e-3.
        @test rel_err(res, ref) < 2e-2
    end

    # All outer paths shared the same inner samples, so their error did not average out (delta -15%, gamma -127% here).
    @testset "vibrato saltando on VG with MC inner integration" begin
        model = variance_gamma(0.2, hyper(100.0, 1.0, 1.0, 0.0), 0.01)
        opt = EuropeanOption(1.0, 80.0)
        ref = pricer(model, ZeroRate(0.01), LewisMethod(400.0, 200000), opt)
        res = pricer(model, ZeroRate(0.01), VibratoSaltando(MonteCarloConfiguration(40_000, 10), VibratoMonteCarloStandard(200)), opt)
        test_value_and_derivatives(res, ref; rtol = 1e-2)
    end

    # The default Analytic grid stopped at ±5 standard deviations and missed part of the score terms in the upper tail (volga -45% here).
    @testset "vibrato saltando VG volga with default Analytic grid" begin
        model = VarianceGammaProcess(hyper(0.2, 1.0, 1.0, 0.0), 0.1, 0.01, Underlying(100.0, 0.01))
        ref = pricer(model, ZeroRate(0.01), LewisMethod(400.0, 200000), eu_opt)
        res = pricer(model, ZeroRate(0.01), VibratoSaltando(MonteCarloConfiguration(100_000, 10; sampling = FinancialMonteCarlo.AntitheticMC()), VibratoMonteCarloAnalytic()), eu_opt)
        test_value_and_derivatives(res, ref; rtol = 3e-2)
    end

    # SobolMC fell into the (then SerialMonteCarloConfig) methods, which matched any sampling, and silently used the rng: on Heston
    # (simulate_with_vol) and for G_T of VG and NIG the results were those of plain MC.
    @testset "Sobol is used on Heston and on the subordinator" begin
        sobol = MonteCarloConfiguration(1_000, 10; sampling = FinancialMonteCarlo.SobolMC())
        plain = MonteCarloConfiguration(1_000, 10)
        model = heston(0.2, 100.0, 0.0)
        (V, X) = VibratoMonteCarlo.simulate_with_vol(model, rfCurve, sobol, 1.0)
        @test X != VibratoMonteCarlo.simulate_with_vol(model, rfCurve, plain, 1.0)[2]
        # The points of FinancialMonteCarlo's Sobol Heston, which floors the variance at 0 instead of eps: 1.3e-8 observed.
        @test X ≈ FinancialMonteCarlo.simulate(model, rfCurve, sobol, 1.0) rtol = 1e-6
        vg = variance_gamma(0.2, 100.0, 0.0)
        @test VibratoMonteCarlo.simulate_subordinator(vg, sobol, 1.0)[1] != VibratoMonteCarlo.simulate_subordinator(vg, plain, 1.0)[1]
    end

    # AntitheticMC drew G_T with rand, as plain MC: no antithetic pairs and no variance reduction on VG and NIG (#29).
    # Now the quantiles at U and 1 - U. Spread over 6 seeds of 10k paths, observed ratios plain / antithetic: 50 on the
    # VG value, 23 on the NIG delta (the NIG value only gains 2).
    @testset "AntitheticMC on the subordinator (#29)" begin
        anti(seed) = MonteCarloConfiguration(10_000, 10; sampling = FinancialMonteCarlo.AntitheticMC(), seed = seed)
        plain(seed) = MonteCarloConfiguration(10_000, 10; seed = seed)
        vg(S0) = VarianceGammaProcess(0.2, -0.15, 0.2, Underlying(S0, 0.0))
        nig(S0) = NormalInverseGaussianProcess(0.2, -0.1, 0.2, Underlying(S0, 0.0))
        for model in (vg(100.0), nig(100.0))
            G_T = VibratoMonteCarlo.simulate_subordinator(model, anti(1), 1.0)[1]
            law = FinancialMonteCarlo.compute_subordinator_law(model, model.κ, 1.0)
            half = div(length(G_T), 2)
            @test VibratoMonteCarlo.Distributions.cdf.(law, G_T[1:half]) .+ VibratoMonteCarlo.Distributions.cdf.(law, G_T[(half+1):end]) ≈ ones(half) rtol = 1e-8
            @test G_T != VibratoMonteCarlo.simulate_subordinator(model, plain(1), 1.0)[1]
        end
        spread(f, mc) = std([f(mc(seed)) for seed = 1:6])
        for method in (mc -> VibratoSaltando(mc, VibratoMonteCarloGaussLegendre()), VibratoMonteCarlo.ConditionalSaltando)
            vg_value(mc) = pricer(vg(100.0), rfCurve, method(mc), eu_opt)
            nig_delta(mc) = pricer(nig(hyper(100.0, 1.0, 1.0, 0.0)), rfCurve, method(mc), eu_opt).epsilon1
            @test spread(vg_value, anti) < spread(vg_value, plain) / 10
            @test spread(nig_delta, anti) < spread(nig_delta, plain) / 5
        end
    end

    # SobolMC and AntitheticMC covered only the diffusion of Kou and Merton: the jumps came from the rng (#30). Now they
    # come from a fixed number of uniforms per path, the Sobol coordinates after the diffusion or the reflections 1 - u,
    # with FinancialMonteCarlo's schemes (#51).
    @testset "Sobol and antithetic on the Kou and Merton jumps (#30)" begin
        kou_30(S0) = KouProcess(0.2, 3.0, 0.3, 10.0, 8.0, Underlying(S0, 0.01))
        merton_30(S0) = MertonProcess(0.2, 3.0, -0.05, 0.1, Underlying(S0, 0.01))
        # The quantiles of the jump laws, against their means on a midpoint grid of the uniforms.
        u_grid = ((1:100_000) .- 0.5) ./ 100_000
        kou_jumps = kou_30(100.0)
        @test mean(FinancialMonteCarlo.jump_quantile.(Ref(kou_jumps), u_grid)) ≈ kou_jumps.p / kou_jumps.λ₊ - (1 - kou_jumps.p) / kou_jumps.λ₋ rtol = 1e-3
        @test mean(FinancialMonteCarlo.jump_quantile.(Ref(merton_30(100.0)), u_grid)) ≈ -0.05 rtol = 1e-3
        # With Sobol the paths no longer depend on the seed, the jumps included.
        sobol(seed) = MonteCarloConfiguration(1_000, 10; sampling = FinancialMonteCarlo.SobolMC(), seed = seed)
        for model in (kou_30(100.0), merton_30(100.0))
            @test VibratoMonteCarlo.simulate_log_with_jump_times(model, rfCurve, sobol(1), 1.0) == VibratoMonteCarlo.simulate_log_with_jump_times(model, rfCurve, sobol(2), 1.0)
        end
        # Kou against Lewis with 64k Sobol points: relative errors 6.3e-4, 3.3e-4, 2.1e-3 on value, delta and gamma
        # (VibratoMonteCarlo's own scheme before #51: 4.6e-4, 8.2e-4, 1.3e-3; with 16k points the gamma of either scheme
        # is 2.2e-3 to 9.4e-3 depending on the point set); the tolerance is 3 times the largest.
        S0 = hyper(100.0, 1.0, 1.0, 0.0)
        res = pricer(kou_30(S0), rfCurve, VibratoSaltando(MonteCarloConfiguration(64_000, 10; sampling = FinancialMonteCarlo.SobolMC()), VibratoMonteCarloGaussLegendre()), eu_opt)
        ref = pricer(kou_30(S0), rfCurve, LewisMethod(400.0, 200000), eu_opt)
        test_value_and_derivatives(res, ref; rtol = 7e-3)
        # Antithetic jumps: spread of the Kou delta over 32 seeds of 10k paths, 3.4e-3 against 6.0e-3 for plain MC
        # (ratio 1.78; 1.68 with VibratoMonteCarlo's own scheme before #51, which 8 seeds cannot tell apart: 1.41 to 2.2);
        # unbiased, the mean error is within 1.1 standard errors.
        anti(seed) = MonteCarloConfiguration(10_000, 10; sampling = FinancialMonteCarlo.AntitheticMC(), seed = seed)
        plain(seed) = MonteCarloConfiguration(10_000, 10; seed = seed)
        kou_delta(mc) = pricer(kou_30(S0), rfCurve, VibratoSaltando(mc, VibratoMonteCarloGaussLegendre()), eu_opt).epsilon1
        @test std([kou_delta(anti(seed)) for seed = 1:32]) < std([kou_delta(plain(seed)) for seed = 1:32]) / 1.4
    end
end
