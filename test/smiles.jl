@isdefined(TEST_UTILS_LOADED) || include("utils.jl")

# All the strikes of a smile share the paths, so on the same seed every entry equals the option priced alone.
@testset "Smiles" begin
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
        ("vibrato, Black-Scholes", opt -> pricer(bs, rfCurve, Vibrato(mc, gl), opt)),
        ("vibrato, FFT density", opt -> pricer(variance_gamma(σ, S0, d), rfCurve, Vibrato(MonteCarloConfiguration(200, 2), analytic), opt)),
        ("VibratoSaltando, finite activity", opt -> pricer(kou(σ, S0, d), rfCurve, VibratoSaltando(mc, gl), opt)),
        ("VibratoSaltando, infinite activity", opt -> pricer(variance_gamma(σ, S0, d), rfCurve, VibratoSaltando(mc, gl), opt)),
        ("PureLRM, Black-Scholes", opt -> pricer(bs, rfCurve, PureLRM(mc, gl), opt)),
        ("PureLRM, FFT density", opt -> pricer(variance_gamma(σ, S0, d), rfCurve, PureLRM(mc, analytic), opt))
    )
    @testset "$name" for (name, price) in cases
        for (smile, opts) in ((eu_smile, eu_opts), (bin_smile, bin_opts))
            @test price(smile) ≈ price.(opts) rtol = 1e-8
        end
    end

    @testset "against the closed form, with derivatives" begin
        @test pricer(bs, rfCurve, PureLRM(mc, gl), bin_smile) ≈ blsbin.(S0, K, r, T, σ, d, isCall) rtol = 1e-9
        delta_smile = ForwardDiff.derivative(x -> pricer(BlackScholesProcess(σ, Underlying(x, d)), rfCurve, PureLRM(mc, gl), eu_smile), S0)
        @test delta_smile ≈ ForwardDiff.derivative(x -> blsprice.(x, K, r, T, σ, d, isCall), S0) rtol = 1e-8
    end
end
