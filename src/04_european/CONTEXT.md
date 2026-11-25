# Stage 04: European pricers

`pricer` methods for European, binary and smile payoffs.

## Inputs

| Load | Why |
|:--|:--|
| `vibrato.jl` | `Vibrato` on Black–Scholes, `ItoProcess`, Heston, other processes (FFT density); `Conditional`; `LRMInterfaceArguments` |
| `pure_lrm.jl` | `PureLRM`: Gaussian density at maturity, or the FFT one (`pure_lrm_fft`, Heston and Lévy) |
| `vibrato_saltando.jl` | `VibratoSaltando` and `ConditionalSaltando`: after the last jump (Kou, Merton), given `G_T` (VG, NIG) |
| `references/european.md` | what each method conditions on, supported derivatives |
| `references/public_api.md` | dispatch on `EuropeanLike`, unsupported payoffs |
| `references/method_limits.md` | accuracy measured against Lewis |
| `references/conventions.md` | `set_seed!` first, drifts from FinancialFFT, discount at the end |

## Process

- A pricer: `set_seed!`, smile to the host, simulate (stage 02), integrate the last step (stage 03), mean over the
  paths, discount with `discount_factor(r, T)` (`exp(-r * T)` for a number, the integral of a curve, #68).
- Drifts from `drift_correction`, never by hand; closed forms from FinancialToolbox.
- Every pricer must stay inferred (`test/quality.jl`).

## Outputs

`pricer(process, rfCurve, method, option)`: a price (a vector for a smile) carrying the derivatives of the inputs.
Checked by `test/vibrato_models.jl`, `test/levy.jl`, `test/conditional.jl`, `test/smiles.jl`, `test/forwarddiff.jl`,
`test/mooncake.jl`, `test/check/accuracy.jl` (opt-in).

## Verify (Kaimon)

Run in a session of this project after a change here (`ex(q = false, ...)`; setup: `using VibratoMonteCarlo, FinancialMonteCarlo, FinancialToolbox, Test`; `bs, und = BlackScholesProcess(0.2), Underlying(100.0, 0.01)`, `m = Vibrato(MonteCarloConfiguration(10_000, 30), VibratoMonteCarloAnalytic(20, -5.0, 5.0))`). Observed 2026-11-21, Julia 1.13.4, FinancialMonteCarlo 0.9.0, FinancialFFT 0.6.3. The AD number types are test dependencies, not loadable from the project environment: derivatives are checked by `run_tests()`.

| Check | Expect |
|:--|:--|
| `@inferred pricer(bs, und, ZeroRate(0.01), m, EuropeanOption(1.0, 80.0))` | `Float64`, 20.78 against `blsprice` 20.98 (Monte Carlo error at 10 000 paths) |
| `@inferred pricer(bs, und, ZeroRate(0.01), m, FinancialMonteCarlo.EuropeanOptionSmile(1.0, [80.0, 100.0, 120.0], [false, true, true]))` | `Vector{Float64}` ≈ [1.199, 7.793, 2.123] against `blsprice` [1.174, 7.886, 2.126] |
| the same pricer on `BlackScholesProcess(0.2f0)`, `Underlying(100.0f0, 0.01f0)`, `ZeroRate(0.01f0)`, `EuropeanOption(1.0f0, 80.0f0)` | `Float32` (Float32 in, Float32 out) |
