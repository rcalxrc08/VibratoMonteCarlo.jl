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

    # TaylorSeries has no eps: the floors of the Heston variance and of the FFT density go through eps_value.
    @testset "eps_value" begin
        for x in (1.0, hyper(1.0, 1.0, 1.0, 0.0), dual(1.0, 1.0), ForwardDiff.Dual(1.0, 1.0), taylor_expand(identity, 1.0, order = 2), taylor_expand(x -> x, [1.0, 2.0]; order = 2)[1])
            @test VibratoMonteCarlo.eps_value(x) == eps(Float64)
        end
        @test VibratoMonteCarlo.eps_value(1.0f0) == eps(Float32)
    end

    # On the same paths the derivatives are exact, differences observed ≤ 3.1e-12 on gamma ≈ 0.02.
    @testset "same as hyper duals" begin
        analytic = VibratoMonteCarloAnalytic()
        am_opt = AmericanOption(T, K, false)
        cases = (
            ("vibrato, Heston", x -> vibrato(heston(σ, x, d), rfCurve, MonteCarloConfiguration(500, 4), eu_opt, vb_mc)),
            ("vibrato, Heston, antithetic", x -> vibrato(heston(σ, x, d), rfCurve, MonteCarloConfiguration(500, 4, FinancialMonteCarlo.AntitheticMC()), eu_opt, vb_mc)),
            ("vibrato_am, Heston", x -> vibrato_am(heston(σ, x, d), rfCurve, MonteCarloConfiguration(500, 4), am_opt, vb_mc)),
            ("vibrato, FFT density", x -> vibrato(kou(σ, x, d), rfCurve, MonteCarloConfiguration(200, 2), eu_opt, analytic)),
            ("pure_lrm, FFT density", x -> pure_lrm(kou(σ, x, d), rfCurve, MonteCarloConfiguration(200, 2), eu_opt, analytic)),
            ("vibrato_am, FFT density", x -> vibrato_am(kou(σ, x, d), rfCurve, MonteCarloConfiguration(200, 4), am_opt, analytic))
        )
        @testset "$name" for (name, price) in cases
            test_value_and_derivatives(price(taylor_expand(identity, S0_v, order = 2)), price(S0_h); rtol = 1e-9)
        end
    end

    @testset "TaylorN against ForwardDiff, Heston" begin
        mc = MonteCarloConfiguration(500, 4)
        price(x) = vibrato(heston(x[2], x[1], d), rfCurve, mc, eu_opt, vb_mc)
        x0 = [S0_v, σ]
        res = price(taylor_expand(x -> x, x0; order = 2))
        @test [res[1][1], res[1][2]] ≈ ForwardDiff.gradient(price, x0) rtol = 1e-10
        @test [2*res[2][1] res[2][2]; res[2][2] 2*res[2][3]] ≈ ForwardDiff.hessian(price, x0) rtol = 1e-9
    end

    # FinancialFFT ≥ 0.2.1 transforms TaylorN and nested ForwardDiff duals. Hessian entries observed ≤ 8e-11 relative to hyper duals.
    @testset "TaylorN and ForwardDiff.hessian, FFT density" begin
        analytic = VibratoMonteCarloAnalytic()
        x0 = [S0_v, σ]
        cases = (("vibrato", x -> vibrato(kou(x[2], x[1], d), rfCurve, MonteCarloConfiguration(200, 2), eu_opt, analytic)), ("pure_lrm", x -> pure_lrm(kou(x[2], x[1], d), rfCurve, MonteCarloConfiguration(200, 2), eu_opt, analytic)))
        @testset "$name" for (name, price) in cases
            # Hyper duals along the directions i and j give the entry (i, j) of the Hessian.
            direction(i, j) = [hyper(x0[k], Float64(k == i), Float64(k == j), 0.0) for k in eachindex(x0)]
            grad = [price(direction(i, i)).epsilon1 for i in eachindex(x0)]
            hess = [price(direction(i, j)).epsilon12 for i in eachindex(x0), j in eachindex(x0)]
            res = price(taylor_expand(x -> x, x0; order = 2))
            @test [res[1][1], res[1][2]] ≈ grad rtol = 1e-10
            @test [2*res[2][1] res[2][2]; res[2][2] 2*res[2][3]] ≈ hess rtol = 1e-9
            @test ForwardDiff.gradient(price, x0) ≈ grad rtol = 1e-10
            @test ForwardDiff.hessian(price, x0) ≈ hess rtol = 1e-9
        end
    end
end
