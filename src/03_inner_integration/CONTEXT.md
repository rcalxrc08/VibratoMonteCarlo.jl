# Stage 03: inner integration of the last step

The vibrato step: the payoff integrated against the density of the last step, derivatives only in the density weight.

## Inputs

| Load | Why |
|:--|:--|
| `lrm.jl` | `log_density`, `v_value`, `v_mod`, `eps_value`, `carries_derivatives`, `SmileStrike`/`map_strikes`, the integrand with its baseline, `lrm_interface_vec!` per config × inner integration (`inner_samples!`, `inner_mean` for `Standard`), `trapezoid_weights`, `singular_points`, `lrm_paths` and the device matrices |
| `lrm_aug.jl` | the FFT density of the last step: `FFTStepLaw`, `fft_step_law`, `value_law`, `lrm_interface_aug!`, `aug_paths` (blocks of `AUG_BLOCK_NODES` on the device), `floor_density` |
| `references/inner_integration.md` | the core idea, the rules, smiles, inference |
| `references/autodiff.md` | `v_value` must strip every level; Mooncake overlays `v_mod` and `carries_derivatives` |
| `references/gpu.md` | the paths × nodes matrices |
| `references/method_limits.md` | bounds and node counts known to fail |

## Process

- Sample points from values only (`v_value(mu) + v_value(sigma) * z`); every derivative flows through `v_mod`.
- Keep the baseline: the value is unchanged, the derivative noise is cut 10-80×.
- `LRMInterfaceArguments`/`IntegrandArguments` keep `@.` broadcasting over the paths only.
- The host loop and the device matrix of a rule must agree (tested to ≤ 4e-11): change both together.
- A smile is integrated strike by strike (`map_strikes`, a `SmileStrike` per strike, no vector of options): on shared
  nodes (fixed-node rules), on the samples drawn once per path (`Standard`), on the density of the path (FFT), one
  column per strike on the device. No array of the integrand per path on the host (#74).
- Uniform grids integrate with `trapezoid_weights` (PatternVectors), never with reductions over the nodes (#69).

## Outputs

`lrm_paths(mu, sigma, args, mcBaseData, vb_mc)` and `aug_paths(law, option, vb_mc, S_end, mcBaseData)`: one undiscounted
value per path (a matrix for a smile). Checked by `test/gauss_legendre.jl`, `test/taylor.jl`, `test/regressions.jl`,
`test/quality.jl` (inference), `test/cuda/runtests.jl`.

## Verify (Kaimon)

Run in a session of this project after a change here (`ex(q = false, ...)`; setup: `using VibratoMonteCarlo, FinancialMonteCarlo, FinancialToolbox, Test`; `bs = BlackScholesProcess(0.2, Underlying(100.0, 0.01))`, `m = Vibrato(MonteCarloConfiguration(10_000, 30), VibratoMonteCarloAnalytic(20, -5.0, 5.0))`). Observed 2026-11-21, Julia 1.13.4, FinancialMonteCarlo 0.9.0, FinancialFFT 0.6.3. The AD number types are test dependencies, not loadable from the project environment: derivatives are checked by `run_tests()`.

| Check | Expect |
|:--|:--|
| `@inferred pricer(bs, ZeroRate(0.01), m, BinaryEuropeanOption(1.0, 80.0))` | `Float64`, 0.8333 against `blsbin` 0.8367 (Monte Carlo error at 10 000 paths) |
| `run_tests()` (gauss_legendre, taylor, regressions, quality) | pass |
