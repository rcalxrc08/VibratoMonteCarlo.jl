# Stage 05: American options

`Vibrato` on an `AmericanPayoff`: Longstaff–Schwartz policy, vibrato cash flows.

## Inputs

| Load | Why |
|:--|:--|
| `vibrato_am.jl` | the frozen-decision pricer (Lévy, Heston: `frozen_laws`, `FrozenFFT`, `FrozenHeston`, `step_payoffs`), the Black–Scholes bridge (`exercise_regions`, `positive_intervals`, `bridge_vibrato`, `bridge_step`, `lrm_segment`) |
| `references/american.md` | the two methods, the half steps, the bridge |
| `references/method_limits.md` | the gamma of the frozen method, the bridge variance in `Nstep` |
| `references/gpu.md` | regressions on the host, the rest on the device |

## Process

- Regressions on values only (`v_value`s); derivatives only through the vibrato of the cash flows.
- Never smooth the exercise decisions.
- Sums over segments with callable structs, not closures in generators (inference on Julia 1.10).

## Outputs

`pricer(process, rfCurve, Vibrato(mc, vb_mc), am_opt)`. Checked by `test/american.jl` (PDE of the Bermudan, calls
without dividends, all backends) and `test/cuda/runtests.jl`.

## Verify (Kaimon)

Run in a session of this project after a change here (`ex(q = false, ...)`; setup: `using VibratoMonteCarlo, FinancialMonteCarlo, FinancialToolbox, Test`; `bs, und = BlackScholesProcess(0.2), Underlying(100.0, 0.01)`, `m = Vibrato(MonteCarloConfiguration(10_000, 30), VibratoMonteCarloAnalytic(20, -5.0, 5.0))`). Observed 2026-11-21, Julia 1.13.4, FinancialMonteCarlo 0.9.0, FinancialFFT 0.6.3. The AD number types are test dependencies, not loadable from the project environment: derivatives are checked by `run_tests()`.

| Check | Expect |
|:--|:--|
| `run_tests()` (american: PDE of the Bermudan, calls without dividends, all backends) | pass |
