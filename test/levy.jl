@isdefined(TEST_UTILS_LOADED) || include("utils.jl")
using ReverseDiff

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

# κ and T change the law of G_T: the draws come from the values, the derivatives from the density of the subordinator.
@testset "Lévy models, derivatives over κ and T" begin
    rfCurve = ZeroRate(0.02)
    var_gamma(κ, S0 = 100.0) = VarianceGammaProcess(0.2, -0.15, κ, Underlying(S0, 0.0))
    nig(κ, S0 = 100.0) = NormalInverseGaussianProcess(0.2, -0.1, κ, Underlying(S0, 0.0))
    hyper_of(x) = hyper(x, 1.0, 1.0, 0.0)
    models = (("Variance Gamma", var_gamma), ("Normal Inverse Gaussian", nig))

    @testset "log_gamma" begin
        log_gamma = VibratoMonteCarlo.log_gamma
        @test log_gamma(5.0) ≈ log(24.0) rtol = 1e-14
        @test log_gamma(0.5) ≈ log(π) / 2 rtol = 1e-14
        @test log_gamma(40.0) ≈ sum(log, 1.0:39.0) rtol = 1e-14
        # Digamma and trigamma at 1, on the number types that have no loggamma too.
        @test ForwardDiff.derivative(log_gamma, 1.0) ≈ -Base.MathConstants.eulergamma rtol = 1e-14
        @test collect(value_and_derivatives(log_gamma(hyper_of(1.0)))) ≈ [0.0, -Base.MathConstants.eulergamma, π^2 / 6] atol = 1e-13
        @test collect(value_and_derivatives(log_gamma(taylor_expand(identity, 1.0, order = 2)))) ≈ [0.0, -Base.MathConstants.eulergamma, π^2 / 6] atol = 1e-13
    end

    @testset "$name, density" for (name, model) in models
        for (κ, t) in ((0.2, 1.0), (0.2, 0.05), (1.5, 3.0)), g in (1e-3, 0.4, 1.0, 5.0)
            law = FinancialMonteCarlo.compute_subordinator_law(model(κ), κ, t)
            @test VibratoMonteCarlo.subordinator_log_density(model(κ), g, t, κ) ≈ VibratoMonteCarlo.Distributions.logpdf(law, g) rtol = 1e-12
        end
    end

    @testset "$name, no derivatives on κ and T" for (name, model) in models
        mc = MonteCarloConfiguration(100, 10)
        @test VibratoMonteCarlo.simulate_subordinator(model(0.2), mc, 1.0)[2] === true
        @test VibratoMonteCarlo.simulate_subordinator(model(0.2), mc, 1)[2] === true
        @test vibrato_saltando(model(0.2), rfCurve, mc, EuropeanOption(1, 100), VibratoMonteCarloGaussLegendre()) == vibrato_saltando(model(0.2), rfCurve, mc, EuropeanOption(1.0, 100.0), VibratoMonteCarloGaussLegendre())
    end

    # Worst error over 8 seeds with 500k paths ≤ 2 standard deviations of the run, the tolerances are 4 of them.
    # For VG the second derivative over κ has a standard deviation of 0.22 against -0.21: it only catches gross errors.
    mc = MonteCarloConfiguration(500_000, 10)
    lewis = LewisMethod(400.0, 200000)
    cases = (("Variance Gamma", var_gamma, (0.007, 0.08, 0.9), (0.007, 0.025, 0.15)), ("Normal Inverse Gaussian", nig, (0.0018, 0.014, 0.11), (0.0018, 0.007, 0.15)))
    @testset "$name against Lewis" for (name, model, tol_κ, tol_T) in cases
        for (model_h, opt, tol) in ((model(hyper_of(0.2)), EuropeanOption(1.0, 100.0), tol_κ), (model(0.2), EuropeanOption(hyper_of(1.0), 100.0), tol_T))
            res = value_and_derivatives(VibratoMonteCarlo.conditional_saltando(model_h, rfCurve, mc, opt))
            ref = value_and_derivatives(pricer(model_h, rfCurve, lewis, opt))
            @test all(abs.(res .- ref) .< tol)
        end
    end

    # G_T is one dimensional: with SobolMode it is the quantile of its law at Sobol points, and the error decreases like 1 / N
    # (16k points: 10-100 times below the spread of plain MC). The points are deterministic, the tolerances are 3 times
    # the errors observed.
    sobol_cases = (("Variance Gamma", var_gamma, ((1e-3, 1e-4, 1e-6), (1e-3, 4e-2, 1.2), (1e-3, 4e-3, 8e-2))), ("Normal Inverse Gaussian", nig, ((3e-4, 6e-5, 2e-6), (3e-4, 8e-3, 0.15), (3e-4, 3e-4, 1e-2))))
    @testset "$name, Sobol" for (name, model, tols) in sobol_cases
        mc_sobol = MonteCarloConfiguration(16_000, 10, FinancialMonteCarlo.SobolMode())
        @test VibratoMonteCarlo.simulate_subordinator(model(0.2), mc_sobol, 1.0)[1] != VibratoMonteCarlo.simulate_subordinator(model(0.2), MonteCarloConfiguration(16_000, 10), 1.0)[1]
        for ((process, opt), tol) in zip(((model(0.2, hyper_of(100.0)), EuropeanOption(1.0, 100.0)), (model(hyper_of(0.2)), EuropeanOption(1.0, 100.0)), (model(0.2), EuropeanOption(hyper_of(1.0), 100.0))), tols)
            res = value_and_derivatives(vibrato_saltando(process, rfCurve, mc_sobol, opt, VibratoMonteCarloGaussLegendre()))
            ref = value_and_derivatives(pricer(process, rfCurve, lewis, opt))
            @test all(abs.(res .- ref) .< tol)
        end
    end

    @testset "$name, AD backends agree" for (name, model) in models
        mc_small = MonteCarloConfiguration(2_000, 10)
        price_κ(κ) = VibratoMonteCarlo.conditional_saltando(model(κ), rfCurve, mc_small, EuropeanOption(1.0, 100.0))
        price_T(T) = VibratoMonteCarlo.conditional_saltando(model(0.2), rfCurve, mc_small, EuropeanOption(T, 100.0))
        for (price, x) in ((price_κ, 0.2), (price_T, 1.0))
            (val, first_der, second_der) = value_and_derivatives(price(hyper_of(x)))
            @test DualNumbers.epsilon(price(dual(x, 1.0))) ≈ first_der rtol = 1e-11
            @test ForwardDiff.derivative(price, x) ≈ first_der rtol = 1e-11
            @test ForwardDiff.derivative(y -> ForwardDiff.derivative(price, y), x) ≈ second_der rtol = 1e-10
            @test collect(value_and_derivatives(price(taylor_expand(identity, x, order = 2)))) ≈ [val, first_der, second_der] rtol = 1e-10
            @test only(ReverseDiff.gradient(y -> price(y[1]), [x])) ≈ first_der rtol = 1e-11
        end
    end

    @testset "$name, smile" for (name, model) in models
        mc_small = MonteCarloConfiguration(500, 4)
        gl = VibratoMonteCarloGaussLegendre()
        K, isCall = [90.0, 100.0, 115.0], [true, false, true]
        smile = ForwardDiff.derivative(κ -> vibrato_saltando(model(κ), rfCurve, mc_small, FinancialMonteCarlo.EuropeanOptionSmile(1.0, K, isCall), gl), 0.2)
        alone = [ForwardDiff.derivative(κ -> vibrato_saltando(model(κ), rfCurve, mc_small, EuropeanOption(1.0, k, c), gl), 0.2) for (k, c) in zip(K, isCall)]
        @test smile ≈ alone rtol = 1e-8
    end
end
