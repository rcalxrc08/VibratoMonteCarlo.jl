# Stage 01: methods and inner integrations

What the user builds: a pricing method (`Vibrato`, `VibratoSaltando`, `PureLRM`) holding a Monte Carlo configuration and
an inner integration of the last step.

## Inputs

| Load | Why |
|:--|:--|
| `utils.jl` | `AbstractVibrato`, the inner integrations (`VibratoMonteCarloStandard`, `Analytic`, `GaussHermite`, `GaussLegendre`) and their `init_lrm_vec`, `rng_type`, `check_argument`/`check_domain`, `CudaConfig`, `to_host`/`on_device`/`normals!`, `path_mean`, `drift_correction`, `analytic_pricer`; rates as numbers or curves (#68): `interval_integral`, `add_drift`, `discount_factor`, `mean_rate`, `drifted_brownian(!)`, `constant_rate` |
| `methods.jl` | the method types, `EuropeanLike`, `check_mode`, `execution_mode`/`with_mode` |
| `references/public_api.md` | exports, which method prices which payoff on which process |
| `references/inner_integration.md` | what each inner integration does, before adding one |
| `references/conventions.md` | input validation (#58), drifts, reuse of the author's packages |

## Process

- A new inner integration is a subtype of `AbstractVibrato` with `init_lrm_vec` here and an `lrm_interface_vec!` (plus
  `device_samples`/`device_option_paths`) in stage 03; it builds its nodes in `rng_type(paths)`.
- A misconfigured method throws `ArgumentError` in its constructor; batched modes are refused there.
- A new method type subtypes `AbstractVibratoMethod` and gets `with_mode`; its `pricer` methods go to stage 04 or 05.

## Outputs

Method values passed to `pricer(process, rfCurve, method, payoff)`. Checked by `test/quality.jl` (validation, modes,
inference, Aqua) and `test/backends.jl` (no ambiguities in `pricer`).

## Verify (Kaimon)

Run in a session of this project after a change here (`ex(q = false, ...)`; setup: `using VibratoMonteCarlo, FinancialMonteCarlo, FinancialToolbox, Test`; `bs, und = BlackScholesProcess(0.2), Underlying(100.0, 0.01)`, `m = Vibrato(MonteCarloConfiguration(10_000, 30), VibratoMonteCarloAnalytic(20, -5.0, 5.0))`). Observed 2026-11-21, Julia 1.13.4, FinancialMonteCarlo 0.9.0, FinancialFFT 0.6.3. The AD number types are test dependencies, not loadable from the project environment: derivatives are checked by `run_tests()`.

| Check | Expect |
|:--|:--|
| `VibratoMonteCarloAnalytic(0)`, `VibratoMonteCarloAnalytic(100, 3.0, -3.0)` | `ArgumentError` from the constructor, before any pricing |
