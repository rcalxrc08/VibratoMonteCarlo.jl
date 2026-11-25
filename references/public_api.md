# Public API: methods, payoffs, dispatch

_The types a user builds and the `pricer` methods they reach._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

## Exports

`src/VibratoMonteCarlo.jl` exports the method types `Vibrato`, `VibratoSaltando`, `PureLRM`, the inner integrations
`VibratoMonteCarloAnalytic`, `VibratoMonteCarloStandard`, `VibratoMonteCarloGaussLegendre`, and re-exports
FinancialMonteCarlo's `pricer`. Not exported: `VibratoMonteCarloGaussHermite`, `VibratoMonteCarlo.Conditional(mc)` (an
analytic-pricer baseline, Black–Scholes only) and `VibratoMonteCarlo.ConditionalSaltando(mc)` (a test oracle, stage 04).

## Methods

The package plugs into FinancialMonteCarlo's interface, `pricer(mcProcess, underlying, rfCurve, method, payoff)` (since 0.7, FinancialMonteCarlo 0.10), as FinancialFFT's
`LewisMethod` does: each estimator is a method type in `src/01_methods/methods.jl`, subtype of
`FinancialMonteCarlo.AbstractMethod` (through `AbstractVibratoMethod`), holding the `MonteCarloConfiguration` (`mc`) and,
where there is one, the inner integration (`vb_mc`). The free functions `vibrato`, `vibrato_saltando`, `vibrato_am`,
`pure_lrm` were removed after 0.2.0.

| Method | Payoffs | Processes | File |
|:--|:--|:--|:--|
| `Vibrato(mc, vb_mc)` | European, binary, smiles | Black–Scholes, `ItoProcess`, Heston, any other (FFT density of the last step) | `src/04_european/vibrato.jl` |
| `VibratoSaltando(mc, vb_mc)` | European, binary, smiles | Kou, Merton (last jump time), VG, NIG (subordinator) | `src/04_european/vibrato_saltando.jl` |
| `PureLRM(mc, vb_mc)` | European, binary, smiles | Gaussian `ItoProcess` (Black–Scholes), Heston and Lévy through the FFT density | `src/04_european/pure_lrm.jl` |
| `Vibrato(mc, vb_mc)` | one `AmericanPayoff` | Black–Scholes (bridge vibrato), Lévy and Heston (frozen decisions) | `src/05_american/vibrato_am.jl` |

`execution_mode`/`with_mode` of a Vibrato method are those of its `mc` (end of `methods.jl`). The batched modes
(`MultiThreading`, `MultiProcess`) are refused when a method is built (`check_mode`: one generator per batch, `host_rng`
errors there).

## Payoffs

- European methods take `EuropeanLike` (an `EuropeanPayoff` or an `AbstractSmile`), so that `Vibrato` on an
  `AmericanPayoff` dispatches to the American method without ambiguities.
- `EuropeanOption`, `BinaryEuropeanOption`, or a smile of them, `EuropeanOptionSmile` or `BinaryEuropeanOptionSmile`
  (FinancialMonteCarlo, not exported), which has **one** maturity and returns a vector of prices.
- Arrays of payoffs, `Position` and `Portfolio` are not supported (FinancialMonteCarlo types those `pricer` methods on
  `MonteCarloConfiguration`).
- Never reimplement payoffs: `payout_f` on a smile calls `FinancialMonteCarlo.payout` through `smile_to_options`, as
  FinancialMonteCarlo's own `pricer` does. That allocates a vector of options at every node; the proper
  `payout(S, smile)` is a TODO of FinancialMonteCarlo, tracked in its CLAUDE.md.
- `Vibrato` on American options takes one FinancialMonteCarlo `AmericanPayoff` (`AmericanOption`,
  `BinaryAmericanOption`); the integrands call `FinancialMonteCarlo.payout`, so any payoff that defines it works.
- `analytic_pricer` (`src/01_methods/utils.jl`) exists only for European and Binary European options.
- The European pricers take a smile to the host first (`to_mode(…, SerialMode())`): its single options are built there
  in every mode.

## Inputs and errors

AD numbers (DualNumbers, HyperDualNumbers, ForwardDiff, ReverseDiff, TaylorSeries) go in as model or market inputs and
the derivatives come out of the price; Mooncake differentiates the plain `Float64` call (`references/autodiff.md`). A
misconfigured method throws `ArgumentError` in its constructor, a value out of its domain `DomainError`
(`references/conventions.md`).
