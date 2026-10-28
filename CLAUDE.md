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
- `backends.jl`: all six extensions give the same BS price, delta and gamma, including nested ForwardDiff (`ForwardDiff.hessian`).
- `vibrato_models.jl`: Kou/BS/Heston/VG × {standard+analytic, standard+MC, antithetic+MC, Sobol+MC}. Compares against Lewis, `blsprice`, or the MC pricer (Heston).
- `forwarddiff.jl`: multi-partial ForwardDiff duals (Merton, BS).
- `taylor.jl`: `Taylor1` inputs with Gauss–Hermite integration.
- `regressions.jl`: one case per fixed bug (Heston antithetic correlation sign, `pure_lrm` default Analytic grid, Kou last-jump-time alignment, shared inner samples, ±5σ default Analytic bounds).
- `levy.jl`: VG and NIG value, delta and gamma against Lewis, for `vibrato_saltando` and the FFT `vibrato`.
- `american.jl`: `vibrato_am`. BS put value, delta and gamma against a Crank–Nicolson PDE of the Bermudan with the same dates (`bermudan_pde`, defined in the file); calls without dividends against the European price (BS and Kou with gamma, Heston value and delta against the MC pricer); all AD backends agree.
- `check/accuracy.jl`: tight `1e-3` checks against Lewis. **Opt-in** (`VIBRATO_TEST_ACCURACY=true`). It currently fails and errors (see Known issues).

`runner_z.jl` is not included: it needs `vibrato_z` from the untracked Zygote experiments.

The suite uses fixed seeds (`set_seed!`), so results are deterministic for a given Julia version. The seeded MersenneTwister stream differs between Julia 1.10 and 1.11+, so calibrate a tolerance on the spread over several seeds (`SerialMode(seed)`), not on one run. If a numerics change trips one, check the error size before loosening it.

Known package bugs are pinned with `@test_broken`, which turns into an error once a bug is fixed. That is the signal to change it to `@test`. There are none at the moment.

Method limits measured against Lewis (they are not bugs, keep them in mind when choosing test settings):
- With `VibratoMonteCarloStandard` (MC inner integration) gammas are noisier than with quadrature: up to ≈ 4% at 100k paths, `Nstep = 300` (values and deltas within 2e-3).
- `VibratoMonteCarloAnalytic` bounds must reach far enough in the upper tail, where a call payoff is large: at ±5σ the VG volga was 3-45% off, at ±7σ (the default) it is within noise. Finer spacing does not help.
- `vibrato_am` on Lévy and Heston freezes the exercise decisions when differentiating (as Pagès, Pironneau, Sall do): value and delta are accurate, but the gamma misses the exercise boundary term and is biased low, below even the European gamma (-23% to -61% on the American put of Longstaff-Schwartz, K = 40, with Black-Scholes). Never add a smoothing of the decisions to fix it.
- `vibrato_am` on Black-Scholes (vibrato over every exercise date) has unbiased gamma, but its variance grows with `Nstep`: the integrand jumps at the exercise boundary by the gap between the exercise value and the path's realized continuation, and the likelihood ratio weights of a bridge of standard deviation `σ √dt / 2` amplify it (gamma standard deviation per 20k-path run ≈ 6% at `Nstep = 10`, 6.5% at 25, 14-19% at 50, K = 40 put).
- The FFT `vibrato` on Lévy processes needs a smooth density of the last step: for VG keep `dt / κ ≳ 2.5` (below 1 the density is singular and the derivatives are garbage), for NIG and Kou gamma degrades with small `dt` (-21% for NIG at `dt / κ = 0.05`, -8% for Kou at 100 steps). For VG and NIG prefer `vibrato_saltando`, which is exact for any `Nstep`.

Known issues in `check/accuracy.jl`:
- `pure_lrm` on non-Ito processes with `Taylor1` inputs throws: `eps(::Type{Taylor1})` in `lrm_aug.jl`, and Interpolations calls `round(::Taylor1, RoundDown)` (TaylorSeries < 0.23 does not define it).

FinancialMonteCarlo 0.3's antithetic Heston simulator (`heston.jl`, `e2 = -(e1 * ρ + ...)`) uses correlation -ρ, the same bug fixed here in `simulate_with_last_vol`. Heston tests therefore use a serial `pricer` reference.

FinancialFFT 0.1.5's TaylorSeries glue calls `get_order`, which was removed in TaylorSeries ≥ 0.23. Tests therefore use hyper-dual Lewis prices as references, never `Taylor1` ones.

## Architecture

Entry point: [src/VibratoMonteCarlo.jl](src/VibratoMonteCarlo.jl). Exports `vibrato`, `vibrato_saltando`, `vibrato_am`, `pure_lrm`, `VibratoMonteCarloAnalytic`, `VibratoMonteCarloStandard`.

### Core idea

1. Simulate the path with ordinary Monte Carlo up to `T - dt` (or up to the last jump time).
2. Price the last step by integrating the payoff against the known transition density, using the **likelihood ratio method (LRM)**. Sample points are built from `v_value(mu) + v_value(sigma) * z`, so the payoff argument carries no derivative. All parameter dependence flows through the density weight.
3. The density weight is `v_mod(log_density)`, where `v_mod(x) = exp(x - v_value(x))`. It has value 1, and its derivatives are the score terms at every order (an automatic Faà di Bruno).
4. The integrand is `(payoff - baseline) * v_mod(log_density) + baseline`, with `baseline` the payoff at the centre of the density (no derivative). The value is unchanged; the score terms have zero mean, so the derivatives lose the baseline times their integration error. This is the vibrato control variate: it cuts the MC-inner gamma noise 10-80×.

All of this sits in [src/metrics/lrm.jl](src/metrics/lrm.jl):
- **`v_value(x)`** strips the derivative part off an AD number. Each AD backend extension **must** implement it, recursively (`v_value(value(x))`), so that nested numbers (`ForwardDiff.hessian`) lose every level: stripping only the outer one leaks the inner derivatives into the payoff argument and gives wrong second derivatives. The base case is `AbstractFloat`.
- **`v_mod(x)`** has a generic `exp`-based fallback. Backends may specialize it for speed or exactness (Dual, Hyper, only over `AbstractFloat` values; nested numbers use the fallback).
- **`lrm_interface_vec!(mu, sigma, args, mcBaseData, vb_mc)`** does the inner integration. It dispatches on the **MC config type × integration method**:
  - With `VibratoMonteCarloStandard` the inner samples are drawn again for every outer path from `mcBaseData.parallelMode.rng` (the `Z` from `init_lrm_vec` is only a buffer). Sharing them across paths turns their error into a bias that does not average out (deltas were 8-15% off).
  - `SerialMonteCarloConfig` × `VibratoMonteCarloStandard`: plain MC average.
  - `SerialAntitheticMonteCarloConfig` × `Standard`: antithetic average.
  - `SerialSobolMonteCarloConfig` × `Standard`: `cld(Nsim_lrm, 2)` points of a 1-D Sobol sequence, randomly shifted modulo 1 per path (Cranley-Patterson), through `norminvcdf`, each paired with its reflection.
  - any config × `VibratoMonteCarloAnalytic`: quadrature on a uniform grid `[x_min, x_max]` of `Npoints` (default 200 on ±7), always in standard deviations of the log increment around its mean. On Gaussian paths that is the standard normal variable. In `lrm_interface_aug!` (FFT density, non-Ito processes) the mean and std come from the FFT density (`DensitySpline`, built by `spline_density`); the grid must stay inside the spline domain (±15 in log-return for `pure_lrm`, ±18 for `vibrato`).
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
- [vibrato_am.jl](src/metrics/vibrato_am.jl), `vibrato_am(process, rfCurve, mc, am_opt::AmericanPayoff, vb_mc)`: exercise at the `Nstep` grid dates, policy from Longstaff–Schwartz (regression on `1, S, S²` of the in-the-money paths, `S` scaled by `S0`, `v_value`s only). Two methods:
  - `BlackScholesProcess`, vibrato over every exercise date. Pass 1 (independent paths, values only, `exercise_regions`) turns the regressions into exercise regions in log price (`positive_intervals`), fixed in price space. Pass 2 simulates `2 Nstep` half steps (same MC configuration, shared rng) so that every exercise date `t_k` is the midpoint of `[s_k, s_{k+1}]`, `s_k = (k - 1/2) dt`, where `log S_{t_k}` is a Brownian bridge (mean of the ends, std `σ √dt / 2`). Backwards along each path, `bridge_vibrato` integrates "exercise if in the region, else continue with the path value" against the bridge density (`lrm_segment`: Gauss–Legendre on each region piece, points from values, derivatives in `v_mod`); the payoff at `T` gets the usual vibrato over the last half step. Since the optimal boundary does not depend on `S0`, a policy fixed in price space gives the `S0` derivatives of every order (the regression error enters them only at second order).
  - Lévy and Heston, frozen decisions. The decision to exercise at `t_k` is taken at `t_{k-1}` and the cash flow is the vibrato expectation `E[payoff(S_{t_k}) | S_{t_{k-1}}]` over one step (`frozen_laws`: FFT `DensitySpline` for Lévy, which needs `VibratoMonteCarloAnalytic`; for Heston a `GaussianStep` with the variance from `simulate_with_vol`). One AD cash flow per path at its exercise date. Other `ItoProcess`es throw.

### Options
- `EuropeanOption`, `BinaryEuropeanOption`, or an `Array` of them that share **one** maturity. `get_unique_maturity` enforces this with `only(unique(...))`.
- `vibrato_am` takes one FinancialMonteCarlo `AmericanPayoff` (`AmericanOption`, `BinaryAmericanOption`); the integrands call `FinancialMonteCarlo.payout`, so any payoff that defines it works.
- Array payoffs return a vector of prices.
- `analytic_pricer` exists only for European and Binary European options.

### Simulation helpers ([src/simulation/](src/simulation/))
These are custom simulators that return extra state the metrics need: the last jump time, the variance path (`simulate_with_vol`; `simulate_with_last_vol` keeps its last column), or the subordinator increments. They allocate matrices with an element type derived from a "zero" built out of all parameters (e.g. `isDualZero`), so AD numbers propagate. Follow this pattern in new simulators.

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
- `vibrato_path.jl`
- `vibrato_zygote*.jl`, `simulate_zygote*.jl`: Zygote reverse-mode experiments with hand-written `rrule`s.

Don't wire them into the module, and don't delete them, unless asked.

## Conventions

- Formatting: JuliaFormatter with `margin = 1200`, so long signatures stay on one line. Also `whitespace_typedefs = true` and `trailing_comma = false`.
- Every public metric calls `FinancialMonteCarlo.set_seed!(mcBaseData)` first, so results are reproducible. Keep this in new metrics.
- The code uses `@muladd` (MuladdMacro), `@views`, and broadcast-heavy `@.` style. Constants come from IrrationalConstants (`log2π`, `sqrt2`, `sqrtπ`).
- Prices are discounted with `exp(-r * T)` at the end of each metric.

## CI and docs

- **GitHub Actions** is the CI:
  - `CI.yml` runs the test suite on every push (Julia 1.10 and latest), with the private registry added through `julia-buildpkg`'s `localregistry` input.
  - `CompatHelper.yml` runs nightly and opens the compat-bump PRs that make up most recent commits.
  - `DeployPage.yml` builds the Franklin/PkgPage site in `page/` (math write-up in `page/index.md`) to `gh-pages`.
- GitLab is no longer used: `.gitlab-ci.yml` (build, test, codecov, benchmark, Documenter site from `docs/make.jl`) and the `origin2` remote are legacy. Push only to `origin` (GitHub).
