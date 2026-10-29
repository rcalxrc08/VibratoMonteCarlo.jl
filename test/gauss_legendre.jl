@isdefined(TEST_UTILS_LOADED) || include("utils.jl")
using ReverseDiff

# pure_lrm on Black-Scholes is a single inner integral, so with VibratoMonteCarloGaussLegendre it matches the closed form
# up to the quadrature error, for every order and number type.
@testset "Gauss-Legendre inner integration" begin
    gl = VibratoMonteCarloGaussLegendre()
    mc = MonteCarloConfiguration(1, 1)
    r, d, σ = 0.02, 0.01, 0.2
    rfCurve = ZeroRate(r)
    bs(σ, S0) = BlackScholesProcess(σ, Underlying(S0, d))

    @testset "pure_lrm against blsprice and blsbin" begin
        # Maximum error observed 1.3e-10 (volga of a call, T = 5): elsewhere ≈ 1e-14.
        for T in (0.01, 1.0, 5.0), K in (40.0, 105.0, 1000.0), isCall in (true, false)
            for (S0, σ_) in ((hyper(100.0, 1.0, 1.0, 0.0), σ), (100.0, hyper(σ, 1.0, 1.0, 0.0)))
                res = value_and_derivatives(pure_lrm(bs(σ_, S0), rfCurve, mc, EuropeanOption(T, K, isCall), gl))
                @test all(isapprox.(res, value_and_derivatives(blsprice(S0, K, r, T, σ_, d, isCall)); atol = 1e-9, rtol = 1e-9))
                res = value_and_derivatives(pure_lrm(bs(σ_, S0), rfCurve, mc, BinaryEuropeanOption(T, K, isCall), gl))
                @test all(isapprox.(res, value_and_derivatives(blsbin(S0, K, r, T, σ_, d, isCall)); atol = 1e-9, rtol = 1e-9))
            end
        end
    end

    @testset "number types" begin
        S0, K, T = 100.0, 105.0, 1.0
        opt = EuropeanOption(T, K)
        price(x) = pure_lrm(bs(x[2], x[1]), rfCurve, mc, opt, gl)
        price_ref(x) = blsprice(x[1], K, r, T, x[2], d)
        @test ForwardDiff.hessian(price, [S0, σ]) ≈ ForwardDiff.hessian(price_ref, [S0, σ]) atol = 1e-9
        @test ReverseDiff.gradient(price, [S0, σ]) ≈ ForwardDiff.gradient(price_ref, [S0, σ]) atol = 1e-11
        @test DualNumbers.epsilon(price([Dual(S0, 1.0), σ])) ≈ DualNumbers.epsilon(price_ref([Dual(S0, 1.0), σ])) atol = 1e-12
        S0_t = taylor_expand(identity, S0, order = 5)
        @test price([S0_t, σ]).coeffs ≈ price_ref([S0_t, σ]).coeffs atol = 1e-12
    end

    @testset "array of options" begin
        S0 = hyper(100.0, 1.0, 1.0, 0.0)
        opts = [EuropeanOption(1.0, 95.0), BinaryEuropeanOption(1.0, 105.0), EuropeanOption(1.0, 120.0, false)]
        res = pure_lrm(bs(σ, S0), rfCurve, mc, opts, gl)
        for (i, opt) in enumerate(opts)
            @test all(isapprox.(value_and_derivatives(res[i]), value_and_derivatives(pure_lrm(bs(σ, S0), rfCurve, mc, opt, gl)); atol = 1e-12))
        end
    end
end
