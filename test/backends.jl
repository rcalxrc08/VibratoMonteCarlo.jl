@isdefined(TEST_UTILS_LOADED) || include("utils.jl")
using ReverseDiff

@testset "AD backends agree" begin
    mc = MonteCarloConfiguration(2_000, 10)
    price(S0) = vibrato(BlackScholesProcess(0.2, Underlying(S0, 0.0)), ZeroRate(0.02), mc, EuropeanOption(1.0, 100.0), VibratoMonteCarloAnalytic(50, -5.0, 5.0))
    (val, delta, gamma) = value_and_derivatives(price(hyper(100.0, 1.0, 1.0, 0.0)))
    @test price(100.0) ≈ val rtol = 1e-12

    res_dual = price(dual(100.0, 1.0))
    @test DualNumbers.value(res_dual) ≈ val rtol = 1e-12
    @test DualNumbers.epsilon(res_dual) ≈ delta rtol = 1e-12
    @test ForwardDiff.derivative(price, 100.0) ≈ delta rtol = 1e-12
    # Nested duals: v_value must strip every level (it stripped only the outer one, gamma was negative).
    @test ForwardDiff.derivative(x -> ForwardDiff.derivative(price, x), 100.0) ≈ gamma rtol = 1e-10
    @test only(ForwardDiff.hessian(x -> price(x[1]), [100.0])) ≈ gamma rtol = 1e-10
    @test only(ReverseDiff.gradient(x -> price(x[1]), [100.0])) ≈ delta rtol = 1e-12
    @test collect(value_and_derivatives(price(taylor_expand(identity, 100.0, order = 2)))) ≈ [val, delta, gamma] rtol = 1e-12
end
