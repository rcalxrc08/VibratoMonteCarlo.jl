@isdefined(TEST_UTILS_LOADED) || include(joinpath(@__DIR__, "..", "utils.jl"))
using Statistics

# Accuracy checks against the Lewis method, opt-in: see runtests.jl.
@testset "accuracy against Lewis method" begin
    S0 = 100.0
    K = 105.0
    r = 0.02
    T = 1.0
    d = 0.0
    σ = 0.2
    lewis = LewisMethod(600.0, 20000)
    mc_lrm = MonteCarloConfiguration(10_000, 300; sampling = FinancialMonteCarlo.AntitheticMC())
    lrm(model, und, rfCurve, mc, opt) = pricer(model, und, rfCurve, PureLRM(mc, VibratoMonteCarloAnalytic(20000, -20.0, 20.0)), opt)
    saltando(model, und, rfCurve, mc, opt) = pricer(model, und, rfCurve, VibratoSaltando(mc, VibratoMonteCarlo.VibratoMonteCarloGaussHermite(200)), opt)
    # Market inputs as functions of the differentiated variable.
    bump_S0(x) = (kou(σ, x, d; λ = 50.0)..., ZeroRate(r))
    bump_r(x) = (kou(σ, S0, d; λ = 50.0)..., ZeroRate(x))

    # PureLRM with the analytic density gives the same number on every seed: a fixed tolerance.
    cases = [(name = "PureLRM, European, S0", inputs = bump_S0, x0 = S0, opt = EuropeanOption(T, K), abs_derivatives = false), (name = "PureLRM, Binary, S0", inputs = bump_S0, x0 = S0, opt = BinaryEuropeanOption(T, K), abs_derivatives = true), (name = "PureLRM, European, r", inputs = bump_r, x0 = r, opt = EuropeanOption(T, K), abs_derivatives = false), (name = "PureLRM, Binary, r", inputs = bump_r, x0 = r, opt = BinaryEuropeanOption(T, K), abs_derivatives = false)]
    @testset "$(c.name)" for c in cases
        ref = pricer(c.inputs(hyper(c.x0, 1.0, 1.0, 0.0))..., lewis, c.opt)
        @testset "$(nameof(typeof(x)))" for x in (taylor_expand(identity, c.x0, order = 5), hyper(c.x0, 1.0, 1.0, 0.0))
            res = lrm(c.inputs(x)..., mc_lrm, c.opt)
            test_value_and_derivatives(res, ref; rtol = 1e-3, abs_derivatives = c.abs_derivatives)
        end
    end

    # VibratoSaltando is a Monte Carlo estimate (#34): with 80k antithetic paths its standard deviation over seeds is
    # 0.62% of the value, 0.29% of delta and 1.8% of gamma (measured on 16 seeds, 2026-11-27), so a fixed 1e-3 fails on
    # noise. The test is on the bias instead: the mean over the seeds within 4 standard errors of Lewis, and the spread
    # below twice the measured one, so that a loss of variance reduction also fails. Taylor1 follows the same paths.
    @testset "vibrato saltando, European, S0" begin
        opt = EuropeanOption(T, K)
        ref = collect(value_and_derivatives(pricer(bump_S0(hyper(S0, 1.0, 1.0, 0.0))..., lewis, opt)))
        mc(seed) = MonteCarloConfiguration(80_000, 100; sampling = FinancialMonteCarlo.AntitheticMC(), seed = seed)
        seeds = 0:15
        estimates = reduce(hcat, [collect(value_and_derivatives(saltando(bump_S0(hyper(S0, 1.0, 1.0, 0.0))..., mc(seed), opt))) for seed in seeds])
        max_spread = [1.2e-2, 6e-3, 3.6e-2]
        @testset "$quantity" for (i, quantity) in enumerate(("value", "delta", "gamma"))
            spread = std(estimates[i, :])
            @test abs(mean(estimates[i, :]) - ref[i]) < 4 * spread / sqrt(length(seeds))
            @test spread / abs(ref[i]) < max_spread[i]
        end
        taylor = collect(value_and_derivatives(saltando(bump_S0(taylor_expand(identity, S0, order = 5))..., mc(first(seeds)), opt)))
        @test taylor ≈ estimates[:, 1] rtol = 1e-10
    end
end
