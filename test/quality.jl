# Aqua.jl (ambiguities, piracy, unbound type parameters, stale dependencies, compat entries) and type stability: every
# pricer is inferred (a filtered generator in a sum, or a vector splatted into a tuple, is not).
using Aqua
@testset "Aqua" begin
    Aqua.test_all(VibratoMonteCarlo)
end

@testset "Inferred pricers" begin
    rfCurve = ZeroRate(0.02)
    und = Underlying(100.0, 0.01)
    opt = EuropeanOption(1.0, 105.0)
    smile = FinancialMonteCarlo.EuropeanOptionSmile(1.0, [95.0, 105.0], [true, false])
    anti = FinancialMonteCarlo.AntitheticMC()
    mc(n, s; sampling = anti) = MonteCarloConfiguration(n, s; sampling = sampling)
    gl, analytic = VibratoMonteCarloGaussLegendre(), VibratoMonteCarloAnalytic()
    cases = [
        ("Vibrato, Black Scholes", BlackScholesProcess(0.2, und), Vibrato(mc(200, 4), gl), opt),
        ("Vibrato, Black Scholes, smile", BlackScholesProcess(0.2, und), Vibrato(mc(200, 4), gl), smile),
        ("Vibrato, Black Scholes, MC inner", BlackScholesProcess(0.2, und), Vibrato(mc(200, 4; sampling = FinancialMonteCarlo.SobolMC()), VibratoMonteCarloStandard(32)), opt),
        ("Vibrato, Heston", heston(0.2, 100.0, 0.01), Vibrato(mc(200, 4), gl), opt),
        ("Vibrato, FFT density", variance_gamma(0.2, 100.0, 0.01), Vibrato(mc(200, 2), analytic), opt),
        ("VibratoSaltando, Kou", kou(0.2, 100.0, 0.01), VibratoSaltando(mc(200, 4), gl), opt),
        ("VibratoSaltando, Merton, smile", merton(0.2, 100.0, 0.01), VibratoSaltando(mc(200, 4), gl), smile),
        ("VibratoSaltando, VG", variance_gamma(0.2, 100.0, 0.01), VibratoSaltando(mc(200, 4), gl), opt),
        ("PureLRM, FFT density", kou(0.2, 100.0, 0.01), PureLRM(mc(200, 2), analytic), opt),
        ("PureLRM, Heston", heston(0.2, 100.0, 0.01), PureLRM(mc(200, 2), analytic), opt),
        ("Vibrato, American", BlackScholesProcess(0.2, und), Vibrato(mc(200, 4), gl), AmericanOption(1.0, 105.0, false)),
    ]
    @testset "$name" for (name, model, method, payoff) in cases
        @test (@inferred pricer(model, rfCurve, method, payoff)) isa Union{Float64, Vector{Float64}}
    end
end

# The mode of a Vibrato method is the one of its configuration (FinancialMonteCarlo 0.6); the batched modes are refused.
@testset "Execution mode" begin
    mc = MonteCarloConfiguration(200, 4; seed = 3)
    gl = VibratoMonteCarloGaussLegendre()
    for method in (Vibrato(mc, gl), VibratoSaltando(mc, gl), PureLRM(mc, gl), VibratoMonteCarlo.Conditional(mc), VibratoMonteCarlo.ConditionalSaltando(mc))
        @test FinancialMonteCarlo.execution_mode(method) == FinancialMonteCarlo.SerialMode()
        moved = FinancialMonteCarlo.with_mode(method, FinancialMonteCarlo.SerialMode())
        @test typeof(moved) == typeof(method) && moved.mc.Nsim == 200 && moved.mc.rng.seed == 3
    end
    @test_throws ArgumentError Vibrato(MonteCarloConfiguration(200, 4; mode = FinancialMonteCarlo.MultiThreading()), gl)
    @test_throws ArgumentError VibratoMonteCarlo.Conditional(MonteCarloConfiguration(200, 4; mode = FinancialMonteCarlo.MultiThreading()))
end

# Invalid inputs throw explicit exceptions, not AssertionError (#58): a misconfigured method fails in its constructor with
# ArgumentError, a maturity out of its domain with DomainError.
@testset "Input validation" begin
    @test_throws ArgumentError VibratoMonteCarloStandard(0)
    @test_throws ArgumentError VibratoMonteCarloAnalytic(0)
    @test_throws ArgumentError VibratoMonteCarloAnalytic(100, 3.0, -3.0)
    @test_throws ArgumentError VibratoMonteCarlo.VibratoMonteCarloGaussHermite(-1)
    @test_throws ArgumentError VibratoMonteCarloGaussLegendre(0)
    model = HestonProcess(0.2, 0.2, 0.0, 0.2, 0.2, 0.2, Underlying(100.0, 0.0))
    for mc in (MonteCarloConfiguration(10, 2), MonteCarloConfiguration(10, 2; sampling = FinancialMonteCarlo.AntitheticMC()), MonteCarloConfiguration(10, 2; sampling = FinancialMonteCarlo.SobolMC()))
        @test_throws DomainError VibratoMonteCarlo.simulate_with_vol(model, ZeroRate(0.01), mc, -1.0)
    end
end

# A Float32 model gives a Float32 price (#45): the inner samples, nodes and grids take the type of the paths.
@testset "Float32" begin
    und = Underlying(100.0f0, 0.0f0)
    rf32 = ZeroRate(0.01f0)
    opt32 = EuropeanOption(1.0f0, 100.0f0)
    bs = BlackScholesProcess(0.2f0, und)
    mc = MonteCarloConfiguration(200, 4)
    anti = MonteCarloConfiguration(200, 4; sampling = FinancialMonteCarlo.AntitheticMC())
    sobol = MonteCarloConfiguration(200, 4; sampling = FinancialMonteCarlo.SobolMC())
    cases = [
        ("Black Scholes, MC inner", bs, Vibrato(mc, VibratoMonteCarloStandard(32))),
        ("Black Scholes, antithetic MC inner", bs, Vibrato(anti, VibratoMonteCarloStandard(32))),
        ("Black Scholes, Sobol MC inner", bs, Vibrato(sobol, VibratoMonteCarloStandard(32))),
        ("Black Scholes, Gauss-Legendre", bs, Vibrato(mc, VibratoMonteCarloGaussLegendre())),
        ("Black Scholes, Gauss-Hermite", bs, Vibrato(mc, VibratoMonteCarlo.VibratoMonteCarloGaussHermite(64))),
        ("Black Scholes, Analytic", bs, Vibrato(mc, VibratoMonteCarloAnalytic())),
        ("Heston, antithetic", HestonProcess(0.2f0, 0.2f0, 0.0f0, 0.2f0, 0.2f0, 0.2f0, und), Vibrato(anti, VibratoMonteCarloGaussLegendre())),
        ("Kou, VibratoSaltando", KouProcess(0.2f0, 3.0f0, 0.4f0, 10.0f0, 5.0f0, und), VibratoSaltando(mc, VibratoMonteCarloGaussLegendre())),
        ("VG, FFT density", VarianceGammaProcess(0.2f0, -0.1f0, 0.2f0, und), Vibrato(mc, VibratoMonteCarloAnalytic())),
        ("PureLRM, Black Scholes", bs, PureLRM(mc, VibratoMonteCarloGaussLegendre())),
    ]
    @testset "$name" for (name, model, method) in cases
        @test pricer(model, rf32, method, opt32) isa Float32
    end
    # the log density of G_T, which weights the derivatives over κ and T of VG and NIG (#70)
    @testset "log_gamma and the subordinator log density" begin
        @test VibratoMonteCarlo.log_gamma(0.5f0) isa Float32
        @test VibratoMonteCarlo.log_gamma(0.5f0) ≈ log(Float32(π)) / 2
        @test VibratoMonteCarlo.log_gamma(20.0f0) isa Float32
        @testset "$(nameof(typeof(model)))" for model in (VarianceGammaProcess(0.2f0, -0.1f0, 0.2f0, und), NormalInverseGaussianProcess(0.2f0, -0.1f0, 0.2f0, und))
            @test @inferred(VibratoMonteCarlo.subordinator_log_density(model, 0.8f0, 1.0f0, model.κ)) isa Float32
        end
    end
end
