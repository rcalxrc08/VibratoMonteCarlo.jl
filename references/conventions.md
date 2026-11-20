# Code conventions

_Rules for every change._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

- Formatting: JuliaFormatter with `margin = 1200`, so long signatures stay on one line. Also
  `whitespace_typedefs = true` and `trailing_comma = false` (`.JuliaFormatter.toml`).
- Every public pricer calls `FinancialMonteCarlo.set_seed!(mcBaseData)` first, so results are reproducible. Keep this
  in new pricers.
- Reuse FinancialMonteCarlo, FinancialFFT and FinancialToolbox before writing a helper: payoffs (`payout`), laws
  (`compute_subordinator_law`), closed forms (`blsprice`, `blsbin`). Write one only when they have nothing generic in the
  number type (`log_gamma`, `subordinator_log_density`: `loggamma` and `Distributions.logpdf` are not).
- Drifts come from FinancialFFT: `drift_correction(process) = FinancialFFT.characteristic_exponent_i(1, process)`, i.e.
  `log E[e^{X_1}]` of the driftless log price, and the risk neutral drift of `log S` is `r - d - drift_correction`.
  Never write `σ² / 2` or `λ * compute_drift` by hand. A process without the two-argument `characteristic_exponent_i`
  (Heston, the log-normal mixtures) throws in the generic `ItoProcess` methods.
- Write `a * b + c` under `@muladd` (or `muladd`) in new code.
- The code uses `@muladd` (MuladdMacro), `@views`, and broadcast-heavy `@.` style. Constants come from
  IrrationalConstants (`log2π`, `sqrt2`, `sqrtπ`).
- Prices are discounted with `exp(-r * T)` at the end of each pricer.
- Input validation (#58, 0.6.0): never `@assert`, which throws a generic `AssertionError` and may be disabled.
  `check_argument` (`ArgumentError`) for a misconfigured method, in its constructor, so it fails before any pricing (and
  before a FinancialCalibrator run, whose loss turns every `DomainError` into `Inf`); `check_domain` (`DomainError`) for
  a value out of its domain (`T > 0` in the Heston simulations). Local copies of FinancialMonteCarlo's helpers
  (`src/01_methods/utils.jl`), not imports of its internals.
- Number types: build arrays in the type of the inputs (`references/autodiff.md`, #45); `eps_value(x)`, not
  `eps(typeof(x))`.
