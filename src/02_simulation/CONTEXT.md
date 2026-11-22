# Stage 02: outer simulation

The paths up to the last step, with the extra state a pricer conditions on.

## Inputs

| Load | Why |
|:--|:--|
| `simulate_with_jump_times.jl` | Kou and Merton: `draw_jumps!` into FinancialMonteCarlo's `JumpList`, `last_jumps`, `simulate_log_with_jump_times` |
| `simulate_with_last_vol.jl` | Heston with the variance path: `simulate_with_vol` (plain, antithetic, Sobol) on FinancialMonteCarlo's Euler kernels (`heston_columns!`, #61), `simulate_with_last_vol` |
| `subordinator.jl` | VG and NIG: `simulate_subordinator` (`G_T` and its weight on `κ`, `T`), `log_gamma`, `subordinator_log_density`, `subordinator_mean` |
| `references/simulation.md` | the conditioning of each simulator, the schemes per sampling, why `G_T` and not the last step |
| `references/gpu.md` | what is drawn on the host and what stays on the device |
| `references/autodiff.md` | element types from a zero of all the parameters |

## Process

- Reuse FinancialMonteCarlo's paths, jumps and subordinator draws; only the extra state is computed here.
- One method per sampling (plain, `AntitheticMC`, `SobolMC`): a missing one silently falls back to the rng.
- Arrays in the type of the parameters (the `isDualZero` pattern), constants in the type of the inputs (#70),
  `eps_value` for epsilons, `@muladd` for `a * b + c`.
- No scalar indexing of device arrays: gathers by linear index, whole-array copies.

## Outputs

Paths and state on the device of the configuration. Checked by `test/regressions.jl` (alignment of the last jump,
antithetic and Sobol on jumps and `G_T`), `test/levy.jl` (`log_gamma`, `subordinator_log_density`),
`test/cuda/runtests.jl`.

## Verify (Kaimon)

Run in a session of this project after a change here (`ex(q = false, ...)`; setup: `using VibratoMonteCarlo, FinancialMonteCarlo, FinancialToolbox, Test`; `bs = BlackScholesProcess(0.2, Underlying(100.0, 0.01))`, `m = Vibrato(MonteCarloConfiguration(10_000, 30), VibratoMonteCarloAnalytic(20, -5.0, 5.0))`). Observed 2026-11-21, Julia 1.13.4, FinancialMonteCarlo 0.9.0, FinancialFFT 0.6.3. The AD number types are test dependencies, not loadable from the project environment: derivatives are checked by `run_tests()`.

| Check | Expect |
|:--|:--|
| `run_tests()` (regressions, levy, CUDA opt-in) | pass; the summary reports Julia 1.13 and `Manifest-v1.13.toml` |
