# European pricers: Vibrato, VibratoSaltando, PureLRM

_Stage 04 reference._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

Every pricer calls `FinancialMonteCarlo.set_seed!(mcBaseData)` first, takes a smile to the host
(`to_mode(…, SerialMode())`), computes `path_mean(lrm_paths(…))` (or `aug_paths`, `subordinator_mean`) and discounts
with `discount_factor(r, T)` at the end. Drifts come from `drift_correction` (`references/conventions.md`). Limits measured against
Lewis: `references/method_limits.md`.

## `src/04_european/vibrato.jl`: `pricer(process, rfCurve, Vibrato(mc, vb_mc), option)`

- `BlackScholesProcess` / `ItoProcess`: Gaussian last step.
- `HestonProcess`: uses `simulate_with_last_vol` to get the last-step variance (stage 02).
- any other process: the FFT density of the last step, `FinancialFFT.density(process, rfCurve, FFTDensityInverter(),
  dt)` wrapped by `fft_step_law`, then `aug_paths` (`src/03_inner_integration/lrm_aug.jl`): `lrm_interface_aug!` path by
  path on the host, a matrix of paths × nodes on the device. The density does not depend on `S_t`: it is evaluated at
  `x - log(S_t)`. This path covers Lévy processes. Interpolations and FFTW are no longer dependencies.
- `VibratoMonteCarlo.Conditional(mc)` (not exported) is an analytic-pricer baseline (BS only).
- `LRMInterfaceArguments(Z, option)` is defined here (`references/inner_integration.md`).

## `src/04_european/vibrato_saltando.jl`: `pricer(process, rfCurve, VibratoSaltando(mc, vb_mc), option)`

- `FiniteActivityProcess` (Kou, Merton): conditions on the **last jump time** (`simulate_log_with_jump_times`, stage
  02) and integrates the Gaussian diffusion from there to `T`. Differentiating w.r.t. `λ` or `T` is unsupported and only
  emits a warning.
- `InfiniteActivityProcess` (VG, NIG): conditions on the whole subordinator path, i.e. on `G_T`
  (`simulate_subordinator`, stage 02), with derivatives over `κ` and `T` through the weight on its law.
- `VibratoMonteCarlo.ConditionalSaltando(mc)` (not exported): same paths as `VibratoSaltando`, but prices what the
  conditioning leaves in closed form (`analytic_pricer`, European and binary only) and differentiates pathwise, with the
  same weight on the law of `G_T` for infinite activity. After the last jump for `FiniteActivityProcess`; for
  `InfiniteActivityProcess` Black-Scholes with volatility `σ √(G_T / T)` and the drift as a dividend yield. It is the
  exact inner-integration limit of `VibratoSaltando`. For infinite activity the price is then a 1-D integral over `G_T`
  (Madan-Carr-Chang for VG), so it is a test oracle more than a pricer: Lewis is the better European pricer.

## `src/04_european/pure_lrm.jl`: `pricer(process, rfCurve, PureLRM(mc, vb_mc), option)`

No path simulation. Pure LRM over the full-maturity density: Gaussian for `ItoProcess` (it needs `σ` and the Lévy drift
correction, so Black-Scholes), the FFT density otherwise and for Heston, whose marginal density at `T` FinancialFFT
computes (`pure_lrm_fft`). The log-normal mixtures have no characteristic exponent in FinancialFFT (FinancialFFT.jl#82),
so `PureLRM`, `Vibrato` and Lewis throw on them. On the device `PureLRM` passes `S0` (FFT) or `mu` (Ito) as a one element
device vector.

## Rates

FinancialMonteCarlo 0.9 curves are isbits and two curves only add or subtract (lazily). Every European method takes a
rate curve and a dividend curve (#68): the drift of the last step and the discount are integrals of the curves, the
Conditional methods price the step with its mean rates, and the FFT density of the last step gets the constant rate
whose drift over the step is the step's (`step_zero_rate`). On the same seed the prices equal those of the constant
rates with the same integrals at `T` (`test/rate_curves.jl`, to 1e-8). Heston's Vibrato (its own Euler scheme) and the
American methods (one frozen law and discount for every step) refuse curves with an `ArgumentError`. Before #68 every
curve case failed with a `MethodError`.
