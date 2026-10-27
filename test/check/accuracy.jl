@isdefined(TEST_UTILS_LOADED) || include(joinpath(@__DIR__, "..", "utils.jl"))

# Tight accuracy checks against the Lewis method, opt-in: see runtests.jl.
@testset "accuracy against Lewis method" begin
    S0 = 100.0
    K = 105.0
    r = 0.02
    T = 1.0
    d = 0.0
    σ = 0.2
    lewis = LewisMethod(600.0, 20000)
    mc_lrm = MonteCarloConfiguration(10_000, 300, FinancialMonteCarlo.AntitheticMC())
    mc_vibrato = MonteCarloConfiguration(80_000, 100, FinancialMonteCarlo.AntitheticMC())
    lrm(model, rfCurve, mc, opt) = pure_lrm(model, rfCurve, mc, opt, VibratoMonteCarloAnalytic(20000, -9.0, 19.0))
    saltando(model, rfCurve, mc, opt) = vibrato_saltando(model, rfCurve, mc, opt, VibratoMonteCarlo.VibratoMonteCarloGaussHermite(200))
    # Market inputs as functions of the differentiated variable.
    bump_S0(x) = (kou(σ, x, d; λ = 50.0), ZeroRate(r))
    bump_r(x) = (kou(σ, S0, d; λ = 50.0), ZeroRate(x))

    cases = [
        (name = "pure_lrm, European, S0", metric = lrm, mc = mc_lrm, inputs = bump_S0, x0 = S0, opt = EuropeanOption(T, K), abs_derivatives = false),
        (name = "pure_lrm, Binary, S0", metric = lrm, mc = mc_lrm, inputs = bump_S0, x0 = S0, opt = BinaryEuropeanOption(T, K), abs_derivatives = true),
        (name = "pure_lrm, European, r", metric = lrm, mc = mc_lrm, inputs = bump_r, x0 = r, opt = EuropeanOption(T, K), abs_derivatives = false),
        (name = "pure_lrm, Binary, r", metric = lrm, mc = mc_lrm, inputs = bump_r, x0 = r, opt = BinaryEuropeanOption(T, K), abs_derivatives = false),
        (name = "vibrato saltando, European, S0", metric = saltando, mc = mc_vibrato, inputs = bump_S0, x0 = S0, opt = EuropeanOption(T, K), abs_derivatives = false),
    ]

    @testset "$(c.name)" for c in cases
        ref = pricer(c.inputs(hyper(c.x0, 1.0, 1.0, 0.0))..., lewis, c.opt)
        @testset "$(nameof(typeof(x)))" for x in (taylor_expand(identity, c.x0, order = 5), hyper(c.x0, 1.0, 1.0, 0.0))
            res = c.metric(c.inputs(x)..., c.mc, c.opt)
            test_value_and_derivatives(res, ref; rtol = 1e-3, abs_derivatives = c.abs_derivatives)
        end
    end
end
