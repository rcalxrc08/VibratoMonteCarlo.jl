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
        mc = MonteCarloConfiguration(20_000, 50, FinancialMonteCarlo.AntitheticMC())
        v, S = VibratoMonteCarlo.simulate_with_last_vol(h, ZeroRate(0.0), mc, 1.0)
        @test cor(v, log.(S[:, end])) < -0.5
    end

    # The Analytic grid was in absolute log-price, the default one truncated at S ≈ 148 (27% too low here).
    @testset "pure_lrm with default Analytic grid on Kou" begin
        model = kou(0.2, 100.0, 0.0)
        ref = pricer(model, rfCurve, LewisMethod(400.0, 200000), eu_opt)
        res = pure_lrm(model, rfCurve, MonteCarloConfiguration(1_000, 10), eu_opt, VibratoMonteCarloAnalytic(2000))
        @test rel_err(res, ref) < 1e-3
    end

    # The last jump time was not aligned with the grid point holding the jump, the drift was counted twice (+7.6% here).
    @testset "vibrato saltando on Kou with coarse grid" begin
        model = kou(0.2, 100.0, 0.0; λ = 50.0)
        ref = pricer(model, rfCurve, LewisMethod(600.0, 20000), eu_opt)
        mc = MonteCarloConfiguration(40_000, 10, FinancialMonteCarlo.AntitheticMC())
        res = vibrato_saltando(model, rfCurve, mc, eu_opt, VibratoMonteCarloAnalytic(200, -6.0, 6.0))
        @test rel_err(res, ref) < 1e-2
    end
end
