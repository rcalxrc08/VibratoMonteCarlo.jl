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
  paths, discount with `exp(-r * T)`.
- Drifts from `drift_correction`, never by hand; closed forms from FinancialToolbox.
- Every pricer must stay inferred (`test/quality.jl`).

## Outputs

`pricer(process, rfCurve, method, option)`: a price (a vector for a smile) carrying the derivatives of the inputs.
Checked by `test/vibrato_models.jl`, `test/levy.jl`, `test/conditional.jl`, `test/smiles.jl`, `test/forwarddiff.jl`,
`test/mooncake.jl`, `test/check/accuracy.jl` (opt-in).
