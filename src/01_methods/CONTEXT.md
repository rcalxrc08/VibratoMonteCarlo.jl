# Stage 01: methods and inner integrations

What the user builds: a pricing method (`Vibrato`, `VibratoSaltando`, `PureLRM`) holding a Monte Carlo configuration and
an inner integration of the last step.

## Inputs

| Load | Why |
|:--|:--|
| `utils.jl` | `AbstractVibrato`, the inner integrations (`VibratoMonteCarloStandard`, `Analytic`, `GaussHermite`, `GaussLegendre`) and their `init_lrm_vec`, `rng_type`, `check_argument`/`check_domain`, `CudaConfig`, `to_host`/`on_device`/`normals!`, `path_mean`, `drift_correction`, `analytic_pricer` |
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
