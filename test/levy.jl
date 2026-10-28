@isdefined(TEST_UTILS_LOADED) || include("utils.jl")

# Value, delta and gamma of infinite activity models against the Lewis method.
@testset "Lévy models against Lewis" begin
    S0 = hyper(100.0, 1.0, 1.0, 0.0)
    rfCurve = ZeroRate(0.02)
    eu_opt = EuropeanOption(1.0, 100.0)
    lewis = LewisMethod(400.0, 200000)
    models = [("Variance Gamma", VarianceGammaProcess(0.2, -0.15, 0.2, Underlying(S0, 0.0))), ("Normal Inverse Gaussian", NormalInverseGaussianProcess(0.2, -0.1, 0.2, Underlying(S0, 0.0)))]

    for (name, model) in models
        ref = pricer(model, rfCurve, lewis, eu_opt)
        @testset "$name" begin
            # Conditioning on the whole subordinator path: exact for any Nstep, errors observed ≤ 9e-4.
            @testset "vibrato saltando" begin
                mc = MonteCarloConfiguration(40_000, 100, FinancialMonteCarlo.AntitheticMC())
                test_value_and_derivatives(vibrato_saltando(model, rfCurve, mc, eu_opt, VibratoMonteCarloAnalytic()), ref; rtol = 5e-3)
            end
            # The FFT density of the last step must be smooth: dt / κ = 2.5 here, errors observed ≤ 2.2e-3.
            @testset "vibrato" begin
                mc = MonteCarloConfiguration(20_000, 2, FinancialMonteCarlo.AntitheticMC())
                test_value_and_derivatives(vibrato(model, rfCurve, mc, eu_opt, VibratoMonteCarloAnalytic()), ref; rtol = 1e-2)
            end
        end
    end
end
