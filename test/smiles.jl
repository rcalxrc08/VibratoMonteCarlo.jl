@isdefined(TEST_UTILS_LOADED) || include("utils.jl")

# All the options of an array or a smile share the paths, so on the same seed every entry equals the option priced alone.
@testset "Arrays of options and smiles" begin
    T, S0, d, r, σ = 1.0, 100.0, 0.01, 0.02, 0.2
    K, isCall = [90.0, 100.0, 115.0], [true, false, true]
    rfCurve = ZeroRate(r)
    eu_smile = FinancialMonteCarlo.EuropeanOptionSmile(T, K, isCall)
    bin_smile = FinancialMonteCarlo.BinaryEuropeanOptionSmile(T, K, isCall)
    eu_opts = EuropeanOption.(T, K, isCall)
    bin_opts = BinaryEuropeanOption.(T, K, isCall)
    gl = VibratoMonteCarloGaussLegendre()
    analytic = VibratoMonteCarloAnalytic()
    mc = MonteCarloConfiguration(500, 4)
    bs = BlackScholesProcess(σ, Underlying(S0, d))

    cases = (
        ("vibrato, Black-Scholes", opt -> vibrato(bs, rfCurve, mc, opt, gl)),
        ("vibrato, FFT density", opt -> vibrato(variance_gamma(σ, S0, d), rfCurve, MonteCarloConfiguration(200, 2), opt, analytic)),
        ("vibrato_saltando, finite activity", opt -> vibrato_saltando(kou(σ, S0, d), rfCurve, mc, opt, gl)),
        ("vibrato_saltando, infinite activity", opt -> vibrato_saltando(variance_gamma(σ, S0, d), rfCurve, mc, opt, gl)),
        ("pure_lrm, Black-Scholes", opt -> pure_lrm(bs, rfCurve, mc, opt, gl)),
        ("pure_lrm, FFT density", opt -> pure_lrm(variance_gamma(σ, S0, d), rfCurve, mc, opt, analytic))
    )
    @testset "$name" for (name, price) in cases
        for (smile, opts) in ((eu_smile, eu_opts), (bin_smile, bin_opts))
            res = price(opts)
            @test res ≈ price.(opts) rtol = 1e-8
            @test price(smile) ≈ res rtol = 1e-14
        end
    end

    @testset "against the closed form, with derivatives" begin
        @test pure_lrm(bs, rfCurve, mc, bin_smile, gl) ≈ blsbin.(S0, K, r, T, σ, d, isCall) rtol = 1e-9
        delta_smile = ForwardDiff.derivative(x -> pure_lrm(BlackScholesProcess(σ, Underlying(x, d)), rfCurve, mc, eu_smile, gl), S0)
        @test delta_smile ≈ ForwardDiff.derivative(x -> blsprice.(x, K, r, T, σ, d, isCall), S0) rtol = 1e-8
    end
end
