# Outer simulation: last jumps, variance paths, subordinator

_Stage 02 reference._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

Files of `src/02_simulation/`: custom simulators that return extra state the pricers need: the last jump time
(`simulate_with_jump_times.jl`), the variance path (`simulate_with_last_vol.jl`: `simulate_with_vol`;
`simulate_with_last_vol` keeps its last column), or the subordinator `G_T` with its weight (`subordinator.jl`:
`simulate_subordinator`). They allocate matrices with an element type derived from a "zero" built out of all parameters
(e.g. `isDualZero`), so AD numbers propagate. Follow this pattern in new simulators. Use `eps_value(x)`, not
`eps(typeof(x))`: TaylorSeries does not define `eps`.

## Last jumps of Kou and Merton (`simulate_with_jump_times.jl`)

- Conditions on the **last jump time** (`simulate_log_with_jump_times`). Jumps are applied at the next grid point, so
  the simulator returns the grid point *before* the step containing the last jump (plus that step's jumps). This keeps
  the conditional step exact and at least one `dt` long; rounding up instead gives `T - T_jump = 0` and NaN derivatives.
- The jumps are FinancialMonteCarlo's (#51): `draw_jumps!` fills a `JumpList` (path, first grid point, size) with its
  `add_jumps!` (plain MC: exponential spacings and `compute_jump_size` from the rng; `AntitheticMC`: the count at a
  uniform `u` and at `1 - u`, then a time `t` and a size uniform `v` per jump, reflected on the path `i + Nsim / 2` as
  the normals are; the sizes by inversion, `jump_quantile`: normal for Merton, double exponential for Kou, continuous at
  `u = 1 - p`), and for `SobolMC` with its `JumpUniforms`, `jump_uniforms` and `add_jumps_from_uniforms!` on the Sobol
  coordinates after the `Nstep` of the diffusion, at the same skip (`FinancialMonteCarlo.sobol_skip!`): each path is one
  point of dimension `Nstep + 1 + 2M` (count at `u[1]`, time and size at `u[2k]`, `u[2k+1]`, `M` the Poisson quantile at
  `1 - sqrt(eps)`, more jumps fall back to the rng).
- `last_jumps` reduces the list in one pass over the jumps (not the paths) to the grid point of the last jump of each
  path and the sum of its jumps; the diffusion at the grid point before it is gathered by linear index, so the paths
  never leave the device.
- Statistically the schemes match the ones of this package before #51 (`JumpUniforms` with sorted times, removed): over
  32 seeds of 10k paths the antithetic Kou delta spread is 1.78 times below plain MC (1.68 before), Sobol on Kou with 64k
  points is within 2.1e-3 of Lewis (1.3e-3 before); with 16k points the Sobol gamma error of either scheme is 2.2e-3 to
  9.4e-3, depending on the point set. With inversion the Kou sizes are differentiable in `p` (pathwise), the rng draws
  are not.

## Heston variance paths (`simulate_with_last_vol.jl`)

`simulate_with_vol` returns `(V, S)`, `V[:, j]` the variance at the time of `S[:, j]`. Serial, Antithetic and Sobol
configs (`simulate_with_vol` has one method each: the generic `MonteCarloConfiguration` one matches any sampling, so a
missing one silently falls back to the rng), `ZeroRate` only, `DomainError` for `T <= 0`. The same code runs on the GPU
(`normals!`). The Euler step is FinancialMonteCarlo's (`heston_log_step`, `heston_variance_step`, through
`heston_columns!`; #61): fused, the antithetic half by the scale `a = -sqrt(dt)` of the increments, the variance floored
at `eps_value` as before (FinancialMonteCarlo floors at `variance_floor`; the floor is outside the kernels, so each
package keeps its own). Same prices as the hand-written step (same seed, to the last digits); measured 2026-11-22 (Julia
1.13.4, Windows laptop; 20000 paths, 50 steps, three interleaved rounds, minimum): ForwardDiff gradient 217-221 ms →
183-196 ms, price within noise (68-80 → 64-75 ms). The antithetic method kept `e2 = e1 ρ + e2_rho sqrt(1 - ρ²)` in a
`Float64` buffer, so a derivative in `ρ` failed with a `MethodError`: fixed by the kernels (test in `vibrato_models.jl`).

## Subordinator of VG and NIG (`subordinator.jl`)

- Given `G_T ~ FinancialMonteCarlo.compute_subordinator_law(process, κ, T)`, `log S_T` is exactly Gaussian with mean
  `log S0 + (r - d - ψ) T + θ G_T` and variance `σ² G_T`, so `Nstep` does not matter. Conditioning on the last step only
  is exact too, but for VG that increment (Gamma, shape `dt / κ`) is often so small that `σ √ΔG` is lost to rounding
  against `log S`.
- `draw_subordinator`: FinancialMonteCarlo's `draw_subordinator!` for plain MC (`rand!` on the host and a kernel on the
  GPU), the quantiles at `U` and `1 - U` for `AntitheticMC`, with the halves ordered as FinancialMonteCarlo's
  `path_sides` (its own antithetic subordinated paths share the subordinator between the halves and reflect the normals,
  which here are integrated out: the halves would be equal), the quantile at 1-D Sobol points for `SobolMC`. The
  quantile of NIG's Inverse Gaussian is a Newton solve, ≈ 7 s for 500k paths against 0.08 s for `rand`, VG's Gamma 1.7 s.
  The path simulators of the subordinated Brownian motion were removed, FinancialMonteCarlo has them.
- `κ` and `T` are supported, unlike `λ`: they change the law of `G_T`, and a draw is not a smooth function of them
  (rejection samplers). `simulate_subordinator` draws `G_T` from their values and returns the weight
  `v_mod(subordinator_log_density(...))`, the LRM step applied to the outer law. Without derivatives on them
  (`carries_derivatives`) the weight is `true`, so no density is computed.
- `subordinator_mean` takes out the sample mean as a baseline, which cuts the noise of the `κ` and `T` derivatives
  3-19×. Its input must be discounted already: a discount outside the baseline multiplies its `T` derivative by the
  score noise.
- `log_gamma` is elementary because DualNumbers, HyperDualNumbers and TaylorSeries have no `loggamma`; the same for
  `subordinator_log_density` (`Distributions.logpdf` is not generic in the number type).
- `draw_subordinator` fills a Float64 buffer: `VibratoSaltando` on VG/NIG returns Float64 prices for Float32 models.
