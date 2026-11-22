# Rate and dividend curves (#68). A European price depends on the rates only through ∫₀ᵀ r and ∫₀ᵀ d, so every method
# priced with curves must give, on the same seed, the price of the constant rates with the same integrals at T: the
# drift of the last step and the discount go through integrals of the curves. Heston's own simulation and the American
# methods use one rate for every step: they refuse curves with an ArgumentError.
using Test, VibratoMonteCarlo, FinancialMonteCarlo

@testset "Rate and dividend curves" begin
    S0, K, T = 100.0, 100.0, 1.0
    r_curve = ZeroRate([0.01, 0.02, 0.04, 0.05], 2.0)
    d_curve = FinancialMonteCarlo.curve_from_nodes([0.0, 0.5, 1.0, 2.0], [0.0, 0.005, 0.01, 0.02])
    r_flat = ZeroRate(FinancialMonteCarlo.integral(r_curve.r, T) / T)
    d_flat = FinancialMonteCarlo.integral(d_curve, T) / T
    mc = MonteCarloConfiguration(2_000, 10; sampling = FinancialMonteCarlo.AntitheticMC(), seed = 1)
    gl, analytic = VibratoMonteCarloGaussLegendre(32), VibratoMonteCarloAnalytic()
    opt = EuropeanOption(T, K)
    bs(d) = BlackScholesProcess(0.2, Underlying(S0, d))
    kou(d) = KouProcess(0.2, 2.0, 0.3, 30.0, 20.0, Underlying(S0, d))
    vg(d) = VarianceGammaProcess(0.2, -0.05, 0.2, Underlying(S0, d))
    heston(d) = HestonProcess(0.2, 0.2, 0.0, 1.5, -0.7, 0.04, Underlying(S0, d))
    # the FFT density of the last step of VG is not smooth at this step (references/method_limits.md): it is left out of
    # the exact comparison, and its grid moves with the drift
    cases = [
        ("Black-Scholes", bs, Vibrato(mc, gl)),
        ("Black-Scholes", bs, PureLRM(mc, gl)),
        ("Black-Scholes", bs, VibratoMonteCarlo.Conditional(mc)),
        ("Kou", kou, Vibrato(mc, analytic)),
        ("Kou", kou, PureLRM(mc, analytic)),
        ("Kou", kou, VibratoSaltando(mc, gl)),
        ("Kou", kou, VibratoMonteCarlo.ConditionalSaltando(mc)),
        ("VarianceGamma", vg, PureLRM(mc, analytic)),
        ("VarianceGamma", vg, VibratoSaltando(mc, gl)),
        ("VarianceGamma", vg, VibratoMonteCarlo.ConditionalSaltando(mc)),
        ("Heston", heston, PureLRM(mc, analytic)),
    ]
    @testset "$name, $(nameof(typeof(method)))" for (name, model, method) in cases
        @test pricer(model(d_curve), r_curve, method, opt) ≈ pricer(model(d_flat), r_flat, method, opt) rtol = 1e-8
    end
    @test_throws ArgumentError pricer(heston(d_flat), r_curve, Vibrato(mc, gl), opt)
    @test_throws ArgumentError pricer(heston(d_curve), r_flat, Vibrato(mc, gl), opt)
    @test_throws ArgumentError pricer(bs(d_flat), r_curve, Vibrato(mc, analytic), AmericanOption(T, K, false))
    @test_throws ArgumentError pricer(kou(d_curve), r_flat, Vibrato(mc, analytic), AmericanOption(T, K, false))
end
