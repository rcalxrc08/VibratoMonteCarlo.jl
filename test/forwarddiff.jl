@isdefined(TEST_UTILS_LOADED) || include("utils.jl")

@testset "ForwardDiff, sensitivities w.r.t. r, S0 and σ" begin
    # Partials are ordered as (r, S0, T, σ), T is not differentiated.
    r = ForwardDiff.Dual(0.02, 1.0, 0.0, 0.0, 0.0)
    S0 = ForwardDiff.Dual(100.0, 0.0, 1.0, 0.0, 0.0)
    σ = ForwardDiff.Dual(0.2, 0.0, 0.0, 0.0, 1.0)
    K = 105.0
    T = 1.0
    d = 0.0
    mc = MonteCarloConfiguration(100_000, 30, FinancialMonteCarlo.AntitheticMC())
    vb_mc = VibratoMonteCarloAnalytic(100, -5.0, 5.0)
    rfCurve = ZeroRate(r)
    eu_opt = EuropeanOption(T, K)
    value_and_partials(x) = [ForwardDiff.value(x); collect(ForwardDiff.partials(x))]
    test_close(res, ref; rtol) = @test all(isapprox.(value_and_partials(res), value_and_partials(ref); rtol = rtol))

    @testset "Merton, vibrato saltando" begin
        model = merton(σ, S0, d)
        res = vibrato_saltando(model, rfCurve, mc, eu_opt, vb_mc)
        test_close(res, pricer(model, rfCurve, mc, eu_opt); rtol = 2e-2)
    end

    @testset "Black Scholes, vibrato" begin
        res = vibrato(BlackScholesProcess(σ, Underlying(S0, d)), rfCurve, mc, eu_opt, vb_mc)
        test_close(res, blsprice(S0, K, r, T, σ, d); rtol = 2e-2)
    end
end
