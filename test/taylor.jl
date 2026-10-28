@isdefined(TEST_UTILS_LOADED) || include("utils.jl")

@testset "TaylorSeries, Gauss-Hermite integration" begin
    S0_v = 100.0
    S0 = taylor_expand(identity, S0_v, order = 5)
    S0_h = hyper(S0_v, 1.0, 1.0, 0.0)
    K = 105.0
    r = 0.02
    T = 1.0
    d = 0.0
    σ = 0.2
    rfCurve = ZeroRate(r)
    eu_opt = EuropeanOption(T, K)
    vb_mc = VibratoMonteCarlo.VibratoMonteCarloGaussHermite()

    @testset "Black Scholes" begin
        # Errors observed over 8 seeds ≤ 1.6e-2 (value), 1.4e-2 (gamma).
        mc = MonteCarloConfiguration(40_000, 30)
        bs_model = BlackScholesProcess(σ, Underlying(S0, d))
        ref = blsprice(S0_h, K, r, T, σ, d)
        @testset "vibrato" begin
            test_value_and_derivatives(vibrato(bs_model, rfCurve, mc, eu_opt, vb_mc), ref; rtol = 3e-2)
        end
        @testset "pure_lrm" begin
            test_value_and_derivatives(pure_lrm(bs_model, rfCurve, mc, eu_opt, vb_mc), ref; rtol = 1e-2)
        end
    end

    @testset "Kou, vibrato saltando" begin
        # Errors observed over 8 seeds ≤ 2.3e-2 (value), 1.3e-2 (delta), 5.1e-2 (gamma).
        mc = MonteCarloConfiguration(40_000, 300, FinancialMonteCarlo.AntitheticMC())
        res = vibrato_saltando(kou(σ, S0, d; λ = 50.0), rfCurve, mc, eu_opt, VibratoMonteCarlo.VibratoMonteCarloGaussHermite(200))
        ref = pricer(kou(σ, S0_h, d; λ = 50.0), rfCurve, LewisMethod(600.0, 20000), eu_opt)
        test_value_and_derivatives(res, ref; rtol = 4e-2, rtol_second = 1e-1)
    end
end
