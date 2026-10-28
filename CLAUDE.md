# CLAUDE.md

Guidance for working in this repository.

## What this is

`VibratoMonteCarlo.jl` is a Julia package that computes option prices **and sensitivities of any order** with Monte Carlo methods that stay smooth (C^∞) under automatic differentiation. You pass AD numbers (DualNumbers, HyperDualNumbers, ForwardDiff, ReverseDiff, TaylorSeries) as model/market inputs, and the derivatives come out of the price.

It builds on the author's own packages:
- **FinancialMonteCarlo.jl**: processes (`BlackScholesProcess`, `KouProcess`, `MertonProcess`, `HestonProcess`, `VarianceGammaProcess`, `NormalInverseGaussianProcess`), payoffs, `MonteCarloConfiguration`, `ZeroRate`, `simulate`, `pricer`, and `set_seed!`.
- **FinancialFFT.jl**: `characteristic_exponent_i`, used to build transition densities for non-Ito processes, and `LewisMethod` reference prices.
- **FinancialToolbox.jl**: `blsprice`/`blsbin`, used for closed-form checks and `analytic_pricer`.

Some dependencies come from a private registry. Add it before `Pkg.instantiate()`:
```julia
import Pkg
Pkg.Registry.add(Pkg.RegistrySpec(url = "https://github.com/rcalxrc08/LocalRegistry.git"))
Pkg.Registry.add("General")
```

## Commands

```bash
julia --project -e 'import Pkg; Pkg.instantiate()'
julia --project -e 'import Pkg; Pkg.test()'                     # runs test/runtests.jl (~12 min on Julia 1.11)
VIBRATO_TEST_ACCURACY=true julia --project -e 'import Pkg; Pkg.test()'   # also runs test/check/accuracy.jl
julia --project=bench -e 'import Pkg; Pkg.develop(path=pwd()); Pkg.instantiate(); include("bench/bench.jl")'
julia -e 'using JuliaFormatter; format(".")'                      # uses .JuliaFormatter.toml
```

`Manifest.toml` files are gitignored, so the first run resolves dependencies fresh.

## Tests

`test/runtests.jl` includes an **explicit list** of files, so new test files must be added there. Each file guards on `TEST_UTILS_LOADED`, so you can also `include` it on its own.

The files:
- `utils.jl`: loads every AD backend, defines model constructors (`kou`, `merton`, `heston`, `variance_gamma`) and `test_value_and_derivatives` (value, 1st and 2nd derivative vs a reference, for `Hyper` and `Taylor1`).
- `backends.jl`: all six extensions give the same BS price, delta and gamma.
- `vibrato_models.jl`: Kou/BS/Heston/VG × {standard+analytic, standard+MC, antithetic+MC, Sobol+MC}. Compares against Lewis, `blsprice`, or the MC pricer (Heston).
- `forwarddiff.jl`: multi-partial ForwardDiff duals (Merton, BS).
- `taylor.jl`: `Taylor1` inputs with Gauss–Hermite integration.
- `regressions.jl`: one case per fixed bug (Heston antithetic correlation sign, `pure_lrm` default Analytic grid, Kou last-jump-time alignment).
- `levy.jl`: VG and NIG value, delta and gamma against Lewis, for `vibrato_saltando` and the FFT `vibrato`.
- `check/accuracy.jl`: tight `1e-3` checks against Lewis. **Opt-in** (`VIBRATO_TEST_ACCURACY=true`). It currently fails and errors (see Known issues).

`runner_z.jl` is not included: it needs `vibrato_z` from the untracked Zygote experiments.

The suite uses fixed seeds (`set_seed!`), so results are deterministic. Tolerances were set from observed errors with roughly 2× margin. If a numerics change trips one, check the error size before loosening it.

Known package bugs are pinned with `@test_broken`, which turns into an error once a bug is fixed. That is the signal to change it to `@test`. There are none at the moment.

Method limits measured against Lewis (they are not bugs, keep them in mind when choosing test settings):
- With `VibratoMonteCarloStandard` (MC inner integration) all outer paths share the same `Nsim_lrm` inner points, so the inner error does not average out: at `Nstep = 300` deltas are ≈ 8% off for every MC config, gammas are unusable. Use `Analytic` or `GaussHermite` when derivatives matter.
- `vibrato_saltando` on Kou/Merton: `T - τ` can be arbitrarily small, so MC inner integration gives very noisy gammas; quadrature is accurate (≤ 0.5%).
- The FFT `vibrato` on Lévy processes needs a smooth density of the last step: for VG keep `dt / κ ≳ 2.5` (below 1 the density is singular and the derivatives are garbage), for NIG and Kou gamma degrades with small `dt` (-21% for NIG at `dt / κ = 0.05`, -8% for Kou at 100 steps). For VG and NIG prefer `vibrato_saltando`, which is exact for any `Nstep`.

Known issues in `check/accuracy.jl`:
- `pure_lrm` on non-Ito processes with `Taylor1` inputs throws: `eps(::Type{Taylor1})` in `lrm_aug.jl`, and Interpolations calls `round(::Taylor1, RoundDown)` (TaylorSeries < 0.23 does not define it).

FinancialMonteCarlo 0.3's antithetic Heston simulator (`heston.jl`, `e2 = -(e1 * ρ + ...)`) uses correlation -ρ, the same bug fixed here in `simulate_with_last_vol`. Heston tests therefore use a serial `pricer` reference.

FinancialFFT 0.1.5's TaylorSeries glue calls `get_order`, which was removed in TaylorSeries ≥ 0.23. Tests therefore use hyper-dual Lewis prices as references, never `Taylor1` ones.

## Architecture

Entry point: [src/VibratoMonteCarlo.jl](src/VibratoMonteCarlo.jl). Exports `vibrato`, `vibrato_saltando`, `pure_lrm`, `VibratoMonteCarloAnalytic`, `VibratoMonteCarloStandard`.

### Core idea

1. Simulate the path with ordinary Monte Carlo up to `T - dt` (or up to the last jump time).
2. Price the last step by integrating the payoff against the known transition density, using the **likelihood ratio method (LRM)**. Sample points are built from `v_value(mu) + v_value(sigma) * z`, so the payoff argument carries no derivative. All parameter dependence flows through the density weight.
3. The density weight is `v_mod(log_density)`, where `v_mod(x) = exp(x - v_value(x))`. It has value 1, and its derivatives are the score terms at every order (an automatic Faà di Bruno).

All of this sits in [src/metrics/lrm.jl](src/metrics/lrm.jl):
- **`v_value(x)`** strips the derivative part off an AD number. Each AD backend extension **must** implement it.
- **`v_mod(x)`** has a generic `exp`-based fallback. Backends may specialize it for speed or exactness (Dual, Hyper).
- **`lrm_interface_vec!(mu, sigma, args, mcBaseData, vb_mc)`** does the inner integration. It dispatches on the **MC config type × integration method**:
  - `SerialMonteCarloConfig` × `VibratoMonteCarloStandard`: plain MC average.
  - `SerialAntitheticMonteCarloConfig` × `Standard`: antithetic average.
  - `SerialSobolMonteCarloConfig` × `Standard`: inner points from a 1-D Sobol sequence through `norminvcdf`, each paired with its reflection (zero mean), built once in `init_lrm_vec`.
  - any config × `VibratoMonteCarloAnalytic`: quadrature on a uniform grid `[x_min, x_max]` of `Npoints`, always in standard deviations of the log increment around its mean. On Gaussian paths that is the standard normal variable. In `lrm_interface_aug!` (FFT density, non-Ito processes) the mean and std come from the FFT density (`DensitySpline`, built by `spline_density`); the grid must stay inside the spline domain (±15 in log-return for `pure_lrm`, ±18 for `vibrato`).
  - `VibratoMonteCarloGaussHermite`: Gauss–Hermite. The type always exists, but it only works once FastGaussQuadrature is loaded (that is what provides `init_lrm_vec`). It is not exported.

### Integration-method types ([src/utils.jl](src/utils.jl))
- `AbstractVibrato` is the supertype.
- Each type implements `init_lrm_vec(vb_mc, mcBaseData)`, which returns the inner sample points or grid.
- `Base.broadcastable` returns `Ref`, so these types are never broadcast over.
- `LRMInterfaceArguments(Z, option)` and `IntegrandArguments(mu, sigma, option)` bundle the vector and option arguments so that `@.` broadcasts only over `mu_jump`/`sigma_jump`. Keep this pattern when adding methods.

### Public metrics (dispatch on process type)
- [vibrato.jl](src/metrics/vibrato.jl), `vibrato(process, rfCurve, mc, option, vb_mc)`:
  - `BlackScholesProcess` / `ItoProcess`: Gaussian last step.
  - `HestonProcess`: uses `simulate_with_last_vol` to get the last-step variance. Serial and Antithetic configs only, `ZeroRate` only.
  - any other process: FFT-derived spline density of the last step (`spline_density` + `lrm_interface_aug!` in `lrm_aug.jl`). This path covers Lévy processes.
  - `conditional(...)` is an analytic-pricer baseline (BS only).
- [vibrato_saltando.jl](src/metrics/vibrato_saltando.jl), `vibrato_saltando(...)`:
  - `FiniteActivityProcess` (Kou, Merton): conditions on the **last jump time** (`simulate_log_with_jump_times`) and integrates the Gaussian diffusion from there to `T`. Jumps are applied at the next grid point, so the simulator returns the grid point *before* the step containing the last jump (plus that step's jumps). This keeps the conditional step exact and at least one `dt` long; rounding up instead gives `T - T_jump = 0` and NaN derivatives. Differentiating w.r.t. `λ` or `T` is unsupported and only emits a warning.
  - `InfiniteActivityProcess` (VG, NIG): conditions on the whole subordinator path, i.e. on `G_T ~ subordinator_law(process, T)`. Given `G_T`, `log S_T` is exactly Gaussian with mean `log S0 + (r - d - ψ) T + θ G_T` and variance `σ² G_T`, so `Nstep` does not matter. Conditioning on the last step only is exact too, but for VG that increment (Gamma, shape `dt / κ`) is often so small that `σ √ΔG` is lost to rounding against `log S`. `simulate_with_sub` is no longer used by the package.
- [pure_lrm.jl](src/metrics/pure_lrm.jl), `pure_lrm(...)`: no path simulation. Pure LRM over the full-maturity density: Gaussian for `ItoProcess`, FFT spline otherwise.

### Options
- `EuropeanOption`, `BinaryEuropeanOption`, or an `Array` of them that share **one** maturity. `get_unique_maturity` enforces this with `only(unique(...))`.
- Array payoffs return a vector of prices.
- `analytic_pricer` exists only for European and Binary European options.

### Simulation helpers ([src/simulation/](src/simulation/))
These are custom simulators that return extra state the metrics need: the last jump time, the last variance, or the subordinator increments. They allocate matrices with an element type derived from a "zero" built out of all parameters (e.g. `isDualZero`), so AD numbers propagate. Follow this pattern in new simulators.

### AD backend glue: package extensions ([ext/](ext/))
- Each AD backend is a Julia package extension (`VibratoMonteCarlo<Pkg>Ext`). The backends are declared under `[weakdeps]`/`[extensions]` in `Project.toml`, which requires Julia ≥ 1.9.
- To support a new number type:
  1. Add `ext/VibratoMonteCarlo<Pkg>Ext.jl` that does `import VibratoMonteCarlo: v_value, v_mod` and adds methods.
  2. Register the package in `[weakdeps]`, `[extensions]` and `[compat]`.
  3. Add it to `test/Project.toml`.
- Extensions cannot define new types in the parent module. `VibratoMonteCarloGaussHermite` therefore lives in `src/utils.jl` and its `lrm_interface_vec!` in `lrm.jl`. Only `init_lrm_vec` (which calls `gausshermite`) is in the FastGaussQuadrature extension.
- ForwardDiff is pulled in transitively (FinancialFFT → DataInterpolations), so its extension is usually active even if the user never loads it.
- `src/deps/zygote_dependencies.jl` is untracked and not wired in (marked `#wrong`).

## Untracked and experimental files — do not treat as package code

`src/VibratoMonteCarlo.jl` does **not** `include` the following files in `src/metrics/`. They are scratch scripts (many are untracked, some are duplicate "copy"/"Copia" files):
- `vibrato_am.jl`: American options via Longstaff–Schwartz on top of `lrm_interface_aug!`. Tracked, WIP.
- `vibrato_path.jl`
- `vibrato_zygote*.jl`, `simulate_zygote*.jl`: Zygote reverse-mode experiments with hand-written `rrule`s.

Don't wire them into the module, and don't delete them, unless asked.

## Conventions

- Formatting: JuliaFormatter with `margin = 1200`, so long signatures stay on one line. Also `whitespace_typedefs = true` and `trailing_comma = false`.
- Every public metric calls `FinancialMonteCarlo.set_seed!(mcBaseData)` first, so results are reproducible. Keep this in new metrics.
- The code uses `@muladd` (MuladdMacro), `@views`, and broadcast-heavy `@.` style. Constants come from IrrationalConstants (`log2π`, `sqrt2`, `sqrtπ`).
- Prices are discounted with `exp(-r * T)` at the end of each metric.

## CI and docs

- **GitLab** (`.gitlab-ci.yml`, Julia 1.10; coverage collects `src/` and `ext/`) is the main CI: build, test with coverage, codecov, benchmark, and a Documenter site (`docs/make.jl`) published to GitLab Pages.
- **GitHub Actions**:
  - `CompatHelper.yml` runs nightly and opens the compat-bump PRs that make up most recent commits.
  - `DeployPage.yml` builds the Franklin/PkgPage site in `page/` (math write-up in `page/index.md`) to `gh-pages`.
- Remotes: `origin` is GitHub and `origin2` is GitLab.
