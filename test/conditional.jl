@isdefined(TEST_UTILS_LOADED) || include("utils.jl")

# conditional_saltando (and conditional for Black-Scholes) prices the last diffusion in closed form: it is the exact inner
# integration limit of vibrato_saltando (vibrato), on the same paths.
@testset "conditional saltando" begin
    rfCurve = ZeroRate(0.02)
    S0 = hyper(100.0, 1.0, 1.0, 0.0)
    mc = MonteCarloConfiguration(10_000, 50)
    vb_mc = VibratoMonteCarlo.VibratoMonteCarloGaussHermite(200)

    for (name, model) in (("Kou", kou(0.2, S0, 0.01)), ("Merton", merton(0.2, S0, 0.01)))
        @testset "$name, European" begin
            opt = EuropeanOption(1.0, 105.0)
            res = VibratoMonteCarlo.conditional_saltando(model, rfCurve, mc, opt)
            # Differences observed over 8 seeds ≤ 1.4e-5 (value), 1.2e-5 (delta), 6e-5 (gamma), the Gauss-Hermite error.
            test_value_and_derivatives(res, vibrato_saltando(model, rfCurve, mc, opt, vb_mc); rtol = 2e-4)
        end
        @testset "$name, binary" begin
            opt = BinaryEuropeanOption(1.0, 105.0)
            res = value_and_derivatives(VibratoMonteCarlo.conditional_saltando(model, rfCurve, mc, opt))
            ref = value_and_derivatives(vibrato_saltando(model, rfCurve, mc, opt, vb_mc))
            # Gauss-Hermite converges slowly on the discontinuous payoff: differences observed over 8 seeds ≤ 1e-3 (value),
            # 1.5e-3 (delta). The gamma is ≈ 1e-4 or less, so it is compared in absolute terms (≤ 6.7e-6 observed).
            @test rel_err(res[1], ref[1]) < 3e-3
            @test rel_err(res[2], ref[2]) < 3e-3
            @test abs(res[3] - ref[3]) < 2e-5
        end
    end

    # Split at the strike, Gauss-Legendre is exact: differences observed ≤ 4.4e-15.
    @testset "exact with Gauss-Legendre" begin
        gl = VibratoMonteCarloGaussLegendre()
        for model in (kou(0.2, S0, 0.01), merton(0.2, S0, 0.01)), opt in (EuropeanOption(1.0, 105.0), BinaryEuropeanOption(1.0, 105.0))
            res = value_and_derivatives(vibrato_saltando(model, rfCurve, mc, opt, gl))
            @test all(isapprox.(res, value_and_derivatives(VibratoMonteCarlo.conditional_saltando(model, rfCurve, mc, opt)); atol = 1e-12))
        end
        bs_model = BlackScholesProcess(0.2, Underlying(S0, 0.01))
        opt = EuropeanOption(1.0, 105.0)
        res = value_and_derivatives(vibrato(bs_model, rfCurve, mc, opt, gl))
        @test all(isapprox.(res, value_and_derivatives(VibratoMonteCarlo.conditional(bs_model, rfCurve, mc, opt)); atol = 1e-12))
    end

    @testset "warns when differentiating over λ" begin
        model = KouProcess(0.2, hyper(5.0, 1.0, 1.0, 0.0), 0.3, 30.0, 20.0, Underlying(100.0, 0.01))
        @test_logs (:warn, r"λ") VibratoMonteCarlo.conditional_saltando(model, rfCurve, MonteCarloConfiguration(100, 10), EuropeanOption(1.0, 105.0))
    end
end
