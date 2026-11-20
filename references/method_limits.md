# Method limits measured against Lewis

_They are not bugs: keep them in mind when choosing test settings or a method._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

- With `VibratoMonteCarloStandard` (MC inner integration) gammas are noisier than with quadrature: up to ≈ 4% at 100k
  paths, `Nstep = 300` (values and deltas within 2e-3).
- `VibratoMonteCarloAnalytic` bounds must reach far enough in the upper tail, where a call payoff is large: at ±5σ the
  VG volga was 3-45% off, at ±7σ (the default) it is within noise. Finer spacing does not help.
- `Vibrato` on American options on Lévy and Heston freezes the exercise decisions when differentiating (as Pagès,
  Pironneau, Sall do): value and delta are accurate, but the gamma misses the exercise boundary term and is biased low,
  below even the European gamma (-23% to -61% on the American put of Longstaff-Schwartz, K = 40, with Black-Scholes).
  Never add a smoothing of the decisions to fix it.
- `Vibrato` on American options on Black-Scholes (vibrato over every exercise date) has unbiased gamma, but its variance
  grows with `Nstep`: the integrand jumps at the exercise boundary by the gap between the exercise value and the path's
  realized continuation, and the likelihood ratio weights of a bridge of standard deviation `σ √dt / 2` amplify it
  (gamma standard deviation per 20k-path run ≈ 6% at `Nstep = 10`, 6.5% at 25, 14-19% at 50, K = 40 put).
- The FFT density path (`aug_paths`: `Vibrato` and `PureLRM` on Lévy processes, `PureLRM` on Heston) integrates with a
  trapezoid that is not split at the strike: binaries converge only O(dx) (`PureLRM` binary call with 200 nodes: 2.4% on
  the value and 11% on gamma for Kou, 2.9% and 46% for Heston; ten times less with 2000), #54. Derivatives that live in
  the tails need wider bounds than ±7: the Heston vol-of-vol derivative of `PureLRM` is 2.9e-3 off on ±7 (2.7e-3 with
  800 nodes), 1.4e-6 on ±10.
- The FFT `Vibrato` on Lévy processes needs a smooth density of the last step: for VG keep `dt / κ ≳ 2.5` (below 1 the
  density is singular and the derivatives are garbage), for NIG and Kou gamma degrades with small `dt` (-21% for NIG at
  `dt / κ = 0.05`, -8% for Kou at 100 steps). For VG and NIG prefer `VibratoSaltando`, which is exact for any `Nstep`.
- `VibratoSaltando` and `ConditionalSaltando` on VG, with `a = T / κ`: at `K = S0 e^{(r - d - ψ) T}` (not the
  at-the-money strike) the n-th `S0` derivative of a call exists only for `a > (n - 1) / 2` and the estimator has finite
  variance only for `a > n - 1` (one order earlier for a binary). This is the smoothness of the true price: the
  conditional gamma scales like `G_T^{-1/2}` and the density of `G_T` like `g^{a - 1}`. `σ`, `θ`, `κ` inherit it through
  `ψ`. Gamma's per-path std is stable in N at `a = 2`, and grows with N below `a = 1` (0.40 to 2.9 from 1e5 to 1e7 paths
  at `a = 0.7`). NIG has no such limit: its density vanishes like `exp(-T² / (2 κ g))` at 0.
- `κ` derivatives through the law of `G_T` are the noisiest (standard deviation at 500k paths, baseline included): VG
  `∂κ` 3%, `∂κ²` 110%; NIG `∂κ` 1.2%, `∂κ²` 12%; `∂T` ≈ 0.1%, `∂T²` ≈ 1.6%.

## Known issues in `test/check/accuracy.jl`

- `VibratoSaltando` on Kou misses the `1e-3` tolerance with 80k paths (0.26% on the value, 0.33% on delta, 2.9% on
  gamma, the same on `Taylor1` and `Hyper`). That is Monte Carlo noise, as the seed spread in `taylor.jl` shows, not a
  defect. The `PureLRM` cases pass.
