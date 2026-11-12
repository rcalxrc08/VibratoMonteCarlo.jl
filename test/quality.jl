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
