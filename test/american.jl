@isdefined(TEST_UTILS_LOADED) || include("utils.jl")
using LinearAlgebra, ReverseDiff

# Reference for vibrato_am under Black-Scholes: the Bermudan option with the same exercise dates (t_1..t_{N-1} and T).
# Crank-Nicolson in log price, two implicit Euler steps after each exercise date. Returns value, delta and gamma at S0.
function bermudan_pde(S0, K, r, d, σ, T, N; iscall = false, nx = 2001, nsub = 40)
    dt = T / N
    L = 10σ * sqrt(T)
    x = range(log(S0) - L, log(S0) + L, length = nx)
    dx = step(x)
    S = exp.(x)
    payoff = iscall ? max.(S .- K, 0.0) : max.(K .- S, 0.0)
    a = σ^2 / 2
    b = r - d - σ^2 / 2
    Lop = Tridiagonal(fill(a / dx^2 - b / (2dx), nx - 1), fill(-2a / dx^2 - r, nx), fill(a / dx^2 + b / (2dx), nx - 1))
    V = copy(payoff)
    hs = dt / nsub
    for e = N:-1:1
        for s = 1:nsub
            θ = s <= 2 ? 1.0 : 0.5
            rhs = V + (1 - θ) * hs * (Lop * V)
            M = Tridiagonal(-θ * hs * Lop.dl, 1 .- θ * hs * Lop.d, -θ * hs * Lop.du)
            M.d[1] = 1
            M.du[1] = 0
            M.d[end] = 1
            M.dl[end] = 0
            τs = (N - e) * dt + s * hs
            rhs[1] = max(blsprice(S[1], K, r, τs, σ, d, iscall), payoff[1])
            rhs[end] = blsprice(S[end], K, r, τs, σ, d, iscall)
            V = M \ rhs
        end
        e > 1 && (V .= max.(V, payoff))
    end
    c = (nx + 1) ÷ 2
    V1 = (V[c+1] - V[c-1]) / (2dx)
    V2 = (V[c+1] - 2V[c] + V[c-1]) / dx^2
    return (V[c], V1 / S0, (V2 - V1) / S0^2)
end

@testset "American options" begin
    K = 40.0
    r = 0.06
    T = 1.0
    σ = 0.2
    S0_v = 40.0
    S0 = hyper(S0_v, 1.0, 1.0, 0.0)
    rfCurve = ZeroRate(r)
    Nstep = 10
    mc = MonteCarloConfiguration(40_000, Nstep, FinancialMonteCarlo.AntitheticMC())
    bs_model = BlackScholesProcess(σ, Underlying(S0, 0.0))

    # Vibrato over every exercise date. Errors observed over 8 seeds ≤ 6.0e-3 (value), 1.1e-2 (delta), 5.4e-2 (gamma);
    # freezing the exercise decisions gives a gamma 40% too low here.
    @testset "Black Scholes put against PDE" begin
        (value, delta, gamma) = bermudan_pde(S0_v, K, r, 0.0, σ, T, Nstep)
        res = vibrato_am(bs_model, rfCurve, mc, AmericanOption(T, K, false), VibratoMonteCarloAnalytic())
        @test rel_err(res.value, value) < 1.5e-2
        @test rel_err(res.epsilon1, delta) < 2.5e-2
        @test rel_err(res.epsilon12, gamma) < 1.2e-1
    end

    # Without dividends early exercise of a call is never optimal: the price is the European one, with all its derivatives.
    @testset "call without dividends is European" begin
        # Errors observed over 8 seeds ≤ 7.6e-3 (value), 9.5e-4 (delta), 7.3e-3 (gamma).
        res = vibrato_am(bs_model, rfCurve, mc, AmericanOption(T, K, true), VibratoMonteCarloAnalytic())
        test_value_and_derivatives(res, blsprice(S0, K, r, T, σ, 0.0); rtol = 1.5e-2, rtol_second = 2e-2)
        # Frozen decisions, FFT density of one step. Errors observed ≤ 1.4e-2 (value), 5.3e-3 (delta), 2.4e-2 (gamma).
        kou_model = kou(σ, S0, 0.0)
        res = vibrato_am(kou_model, rfCurve, MonteCarloConfiguration(20_000, 10, FinancialMonteCarlo.AntitheticMC()), AmericanOption(T, K, true), VibratoMonteCarloAnalytic())
        test_value_and_derivatives(res, pricer(kou_model, rfCurve, LewisMethod(400.0, 200000), EuropeanOption(T, K)); rtol = 3e-2, rtol_second = 5e-2)
        # Frozen decisions, Gaussian step with the simulated variance. The MC reference has no reliable second derivative.
        # Errors observed ≤ 9.8e-3 (value), 4.1e-3 (delta).
        heston_model = heston(σ, S0, 0.0)
        res = vibrato_am(heston_model, rfCurve, MonteCarloConfiguration(40_000, 25, FinancialMonteCarlo.AntitheticMC()), AmericanOption(T, K, true), VibratoMonteCarloAnalytic())
        heston_ref = pricer(heston_model, rfCurve, MonteCarloConfiguration(400_000, 25), EuropeanOption(T, K))
        @test rel_err(res.value, heston_ref.value) < 2e-2
        @test rel_err(res.epsilon1, heston_ref.epsilon1) < 1e-2
    end

    @testset "AD backends agree" begin
        mc_small = MonteCarloConfiguration(2_000, 10)
        price(x) = vibrato_am(BlackScholesProcess(σ, Underlying(x, 0.0)), rfCurve, mc_small, AmericanOption(T, K, false), VibratoMonteCarloAnalytic())
        (val, delta, gamma) = value_and_derivatives(price(S0))
        @test price(S0_v) ≈ val rtol = 1e-12
        @test DualNumbers.epsilon(price(dual(S0_v, 1.0))) ≈ delta rtol = 1e-12
        @test ForwardDiff.derivative(price, S0_v) ≈ delta rtol = 1e-12
        @test only(ForwardDiff.hessian(x -> price(x[1]), [S0_v])) ≈ gamma rtol = 1e-10
        @test only(ReverseDiff.gradient(x -> price(x[1]), [S0_v])) ≈ delta rtol = 1e-12
        @test collect(value_and_derivatives(price(taylor_expand(identity, S0_v, order = 2)))) ≈ [val, delta, gamma] rtol = 1e-12
    end
end
