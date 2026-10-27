using Test, VibratoMonteCarlo, FinancialMonteCarlo, FinancialFFT, FinancialToolbox
# Loading the AD backends triggers the corresponding Requires.jl glue in VibratoMonteCarlo.
using DualNumbers, HyperDualNumbers, ForwardDiff, TaylorSeries, FastGaussQuadrature

const TEST_UTILS_LOADED = true

# Model parameters shared by all tests, only the differentiated inputs vary.
kou(σ, S0, d; λ = 5.0) = KouProcess(σ, λ, 0.3, 30.0, 20.0, Underlying(S0, d))
merton(σ, S0, d) = MertonProcess(σ, 5.0, 0.03, 0.02, Underlying(S0, d))
heston(σ, S0, d) = HestonProcess(σ, 0.2, 0.0, 0.2, 0.2, 0.2, Underlying(S0, d))
variance_gamma(σ, S0, d) = VarianceGammaProcess(σ, 0.1, 0.01, Underlying(S0, d))

rel_err(x, ref) = abs(x - ref) / abs(ref)

# Value, first and second derivative w.r.t. the differentiated input.
value_and_derivatives(x::Hyper) = (x.value, x.epsilon1, x.epsilon12)
value_and_derivatives(x::Taylor1) = (x[0], x[1], 2 * x[2])

"""
    test_value_and_derivatives(result, reference; rtol, rtol_second = rtol, abs_derivatives = false)

Compare value, first and second derivative of `result` against `reference`.
The second derivative is checked against `rtol_second`.
The value is compared in relative terms, the derivatives in absolute terms if `abs_derivatives`.
"""
function test_value_and_derivatives(result, reference; rtol, rtol_second = rtol, abs_derivatives = false)
    res = value_and_derivatives(result)
    ref = value_and_derivatives(reference)
    @test rel_err(res[1], ref[1]) < rtol
    for (i, tol) in ((2, rtol), (3, rtol_second))
        err = abs_derivatives ? abs(res[i] - ref[i]) : rel_err(res[i], ref[i])
        @test err < tol
    end
end
