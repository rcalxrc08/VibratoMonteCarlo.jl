@isdefined(TEST_UTILS_LOADED) || include("utils.jl")
using Mooncake

# The first call compiles Mooncake's rules, which takes minutes: few paths, and only two reverse and one forward case.
@testset "Mooncake" begin
    x0 = [100.0, 0.2, 0.02, 0.01] # S0, σ, r, d
    gradient(f, cache) = copy(Mooncake.value_and_gradient!!(cache, f, x0)[2][2])
    forward_config = Mooncake.Config(enable_nfwd = false)

    @testset "reverse, Gaussian last step" begin
        # A binary has no pathwise derivative: all of it comes from the score terms, which need the v_value and v_mod rules.
        f(x) = vibrato(BlackScholesProcess(x[2], Underlying(x[1], x[4])), ZeroRate(x[3]), MonteCarloConfiguration(200, 5), BinaryEuropeanOption(1.0, 105.0), VibratoMonteCarloGaussLegendre())
        @test gradient(f, Mooncake.prepare_gradient_cache(f, x0)) ≈ ForwardDiff.gradient(f, x0) rtol = 1e-12
    end

    @testset "reverse, law of the subordinator" begin
        # Mooncake sees floats: κ and T need their score terms even though carries_derivatives(::AbstractFloat) is false for the other backends.
        y0 = [x0; 0.2; 1.0] # κ, T
        f(y) = vibrato_saltando(VarianceGammaProcess(y[2], -0.15, y[5], Underlying(y[1], y[4])), ZeroRate(y[3]), MonteCarloConfiguration(200, 5), EuropeanOption(y[6], 105.0), VibratoMonteCarloGaussLegendre())
        @test copy(Mooncake.value_and_gradient!!(Mooncake.prepare_gradient_cache(f, y0), f, y0)[2][2]) ≈ ForwardDiff.gradient(f, y0) rtol = 1e-12
    end

    @testset "forward, FFT density" begin
        f(x) = vibrato(variance_gamma(x[2], x[1], x[4]), ZeroRate(x[3]), MonteCarloConfiguration(200, 2), EuropeanOption(1.0, 105.0), VibratoMonteCarloAnalytic())
        reference = [HyperDualNumbers.eps1(f(hyper.(x0, Float64.(eachindex(x0) .== i), Float64.(eachindex(x0) .== i), 0.0))) for i in eachindex(x0)]
        @test gradient(f, Mooncake.prepare_derivative_cache(f, x0; config = forward_config)) ≈ reference rtol = 1e-12
        @test ForwardDiff.gradient(f, x0) ≈ reference rtol = 1e-12 # FinancialFFT < 0.2 conjugated the partials of fft!
    end
end
