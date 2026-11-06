# CLAUDE.md

Guidance for working in this repository.

## What this is

`VibratoMonteCarlo.jl` is a Julia package that computes option prices **and sensitivities of any order** with Monte Carlo methods that stay smooth (C^∞) under automatic differentiation. You pass AD numbers (DualNumbers, HyperDualNumbers, ForwardDiff, ReverseDiff, TaylorSeries) as model/market inputs, and the derivatives come out of the price. Mooncake is also supported: it differentiates the plain `Float64` call.

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
julia --project=test/cuda -e 'import Pkg; Pkg.instantiate(); include("test/cuda/runtests.jl")'   # GPU tests, needs an NVIDIA GPU (~2 min on an MX150)
```

`Manifest.toml` files are gitignored, so the first run resolves dependencies fresh.

## Tests

`test/runtests.jl` includes an **explicit list** of files, so new test files must be added there. Each file guards on `TEST_UTILS_LOADED`, so you can also `include` it on its own.

The files:
- `utils.jl`: loads every AD backend, defines model constructors (`kou`, `merton`, `heston`, `variance_gamma`) and `test_value_and_derivatives` (value, 1st and 2nd derivative vs a reference, for `Hyper` and `Taylor1`).
- `backends.jl`: all five AD extensions give the same BS price, delta and gamma, including nested ForwardDiff (`ForwardDiff.hessian`). No method ambiguities in `pricer` (`Test.detect_ambiguities`).
- `vibrato_models.jl`: Kou/BS/Heston/VG × {standard+analytic, standard+MC, antithetic+MC, Sobol+MC}. Compares against Lewis, `blsprice`, or the MC pricer (Heston).
- `forwarddiff.jl`: multi-partial ForwardDiff duals (Merton, BS).
- `taylor.jl`: `Taylor1` inputs with Gauss–Hermite integration, against references; on Heston and the FFT density against hyper duals on the same paths; `eps_value`; `TaylorN` on Heston against ForwardDiff; `TaylorN` and `ForwardDiff.hessian` through the FFT density (Kou) against hyper duals.
- `regressions.jl`: one case per fixed bug (Heston antithetic correlation sign, Sobol ignored on Heston and on `G_T`, `PureLRM` default Analytic grid, Kou last-jump-time alignment, shared inner samples, ±5σ default Analytic bounds).
- `levy.jl`: VG and NIG value, delta and gamma against Lewis, for `VibratoSaltando` and the FFT `Vibrato`. Derivatives over `κ` and `T`: `log_gamma`, `subordinator_log_density` against `logpdf`, `ConditionalSaltando` against Lewis, all AD backends agree, a smile against its options alone, and no weight without derivatives.
- `american.jl`: `Vibrato` on American options. BS put value, delta and gamma against a Crank–Nicolson PDE of the Bermudan with the same dates (`bermudan_pde`, defined in the file); calls without dividends against the European price (BS and Kou with gamma, Heston value and delta against the MC pricer); all AD backends agree.
- `conditional.jl`: `ConditionalSaltando` (Kou, Merton, VG, NIG; European and binary) against `VibratoSaltando` on the same paths (Gauss–Hermite within its quadrature error, Gauss–Legendre to 1e-12, on VG and NIG also over `κ` and `T`), `Conditional` against `Vibrato` on BS, and the warning when differentiating over `λ`.
- `gauss_legendre.jl`: `VibratoMonteCarloGaussLegendre` through `PureLRM` on BS against `blsprice`/`blsbin` to 1e-9 (calls, puts, binaries, short and long maturities, strikes outside the bounds), with ForwardDiff Hessians, ReverseDiff, DualNumbers, `Taylor1` and a smile.
- `smiles.jl`: European and binary smiles through `Vibrato`, `VibratoSaltando` and `PureLRM`: every entry against the option priced alone on the same seed, a binary smile against `blsbin`.
- `mooncake.jl`: reverse gradients of a BS binary and of `VibratoSaltando` on VG (with `κ` and `T`) against ForwardDiff, forward gradient of the FFT `Vibrato` on VG against hyper duals (ForwardDiff against them too). Compiling the rules takes minutes, so it has few paths and few cases.
- `check/accuracy.jl`: tight `1e-3` checks against Lewis. **Opt-in** (`VIBRATO_TEST_ACCURACY=true`). It currently fails (see Known issues).
- `examples.jl`: README.md is up to date with the examples (`etc/readme.jl` copies `examples/<name>.jl` between `<!-- BEGIN examples/<name>.jl -->` and `<!-- END ... -->`: GitHub cannot include a file, so after changing an example run `julia etc/readme.jl`), and `examples/hello_world.jl` runs and matches `blsprice`.
- `cuda/runtests.jl`: GPU tests, **opt-in** and in their own project (`test/cuda/Project.toml`, with `[sources]` for the package, Julia ≥ 1.11), so CI neither downloads nor compiles CUDA. With `CUDA.allowscalar(false)`: the device inner integration against the host one on the same `mu` (≤ 4e-11), the GPU Heston paths against FinancialMonteCarlo's GPU Heston (1.4e-8, the variance floors differ), VG/NIG on the same `G_T` as the host (1e-15), the CURAND estimators against references (tolerances twice the worst of 6 seeds), nested ForwardDiff on the device against hyper duals, and the host parts.

`runner_z.jl` is not included: it needs `vibrato_z` from the untracked Zygote experiments.

The suite uses fixed seeds (`set_seed!`), so results are deterministic for a given Julia version. The seeded MersenneTwister stream differs between Julia 1.10 and 1.11+, so calibrate a tolerance on the spread over several seeds (`SerialMode(seed)`), not on one run. If a numerics change trips one, check the error size before loosening it.

Known package bugs are pinned with `@test_broken`, which turns into an error once a bug is fixed. That is the signal to change it to `@test`. There are none at the moment.

Method limits measured against Lewis (they are not bugs, keep them in mind when choosing test settings):
- With `VibratoMonteCarloStandard` (MC inner integration) gammas are noisier than with quadrature: up to ≈ 4% at 100k paths, `Nstep = 300` (values and deltas within 2e-3).
- `VibratoMonteCarloAnalytic` bounds must reach far enough in the upper tail, where a call payoff is large: at ±5σ the VG volga was 3-45% off, at ±7σ (the default) it is within noise. Finer spacing does not help.
- `Vibrato` on American options on Lévy and Heston freezes the exercise decisions when differentiating (as Pagès, Pironneau, Sall do): value and delta are accurate, but the gamma misses the exercise boundary term and is biased low, below even the European gamma (-23% to -61% on the American put of Longstaff-Schwartz, K = 40, with Black-Scholes). Never add a smoothing of the decisions to fix it.
- `Vibrato` on American options on Black-Scholes (vibrato over every exercise date) has unbiased gamma, but its variance grows with `Nstep`: the integrand jumps at the exercise boundary by the gap between the exercise value and the path's realized continuation, and the likelihood ratio weights of a bridge of standard deviation `σ √dt / 2` amplify it (gamma standard deviation per 20k-path run ≈ 6% at `Nstep = 10`, 6.5% at 25, 14-19% at 50, K = 40 put).
- The FFT `Vibrato` on Lévy processes needs a smooth density of the last step: for VG keep `dt / κ ≳ 2.5` (below 1 the density is singular and the derivatives are garbage), for NIG and Kou gamma degrades with small `dt` (-21% for NIG at `dt / κ = 0.05`, -8% for Kou at 100 steps). For VG and NIG prefer `VibratoSaltando`, which is exact for any `Nstep`.
- `VibratoSaltando` and `ConditionalSaltando` on VG, with `a = T / κ`: at `K = S0 e^{(r - d - ψ) T}` (not the at-the-money strike) the n-th `S0` derivative of a call exists only for `a > (n - 1) / 2` and the estimator has finite variance only for `a > n - 1` (one order earlier for a binary). This is the smoothness of the true price: the conditional gamma scales like `G_T^{-1/2}` and the density of `G_T` like `g^{a - 1}`. `σ`, `θ`, `κ` inherit it through `ψ`. Gamma's per-path std is stable in N at `a = 2`, and grows with N below `a = 1` (0.40 to 2.9 from 1e5 to 1e7 paths at `a = 0.7`). NIG has no such limit: its density vanishes like `exp(-T² / (2 κ g))` at 0.
- `κ` derivatives through the law of `G_T` are the noisiest (standard deviation at 500k paths, baseline included): VG `∂κ` 3%, `∂κ²` 110%; NIG `∂κ` 1.2%, `∂κ²` 12%; `∂T` ≈ 0.1%, `∂T²` ≈ 1.6%.

Known issues in `check/accuracy.jl`:
- `VibratoSaltando` on Kou misses the `1e-3` tolerance with 80k paths (0.26% on the value, 0.33% on delta, 2.9% on gamma, the same on `Taylor1` and `Hyper`). That is Monte Carlo noise, as the seed spread in `taylor.jl` shows, not a defect. The `PureLRM` cases pass.

Known gaps of the number types:
- The FFT density path (`Vibrato` and `PureLRM` on non-Ito processes, European and American) takes `TaylorN` and nested ForwardDiff duals (`ForwardDiff.hessian`) because `[compat]` requires FinancialFFT 0.3 (the support came in 0.2.1: with 0.2.0 `plan_fft!` had no method for them). Hessians match hyper duals to ≈ 1e-10 relative on VG, NIG, Kou and Merton.
- `Vibrato` through the FFT density (European and American) does not take derivatives over `κ` on VG and NIG, for any number type: the path simulation builds `Gamma` with AD parameters. `PureLRM` does, and so does `VibratoSaltando`.

## Architecture

Entry point: [src/VibratoMonteCarlo.jl](src/VibratoMonteCarlo.jl). Exports the method types `Vibrato`, `VibratoSaltando`, `PureLRM`, the inner integrations `VibratoMonteCarloAnalytic`, `VibratoMonteCarloStandard`, `VibratoMonteCarloGaussLegendre`, and re-exports FinancialMonteCarlo's `pricer`.

### Core idea

1. Simulate the path with ordinary Monte Carlo up to `T - dt` (or up to the last jump time).
2. Price the last step by integrating the payoff against the known transition density, using the **likelihood ratio method (LRM)**. Sample points are built from `v_value(mu) + v_value(sigma) * z`, so the payoff argument carries no derivative. All parameter dependence flows through the density weight.
3. The density weight is `v_mod(log_density)`, where `v_mod(x) = exp(x - v_value(x))`. It has value 1, and its derivatives are the score terms at every order (an automatic Faà di Bruno).
4. The integrand is `(payoff - baseline) * v_mod(log_density) + baseline`, with `baseline` the payoff at the centre of the density (no derivative). The value is unchanged; the score terms have zero mean, so the derivatives lose the baseline times their integration error. This is the vibrato control variate: it cuts the MC-inner gamma noise 10-80×.

All of this sits in [src/metrics/lrm.jl](src/metrics/lrm.jl):
- **`v_value(x)`** strips the derivative part off an AD number. Each AD backend extension **must** implement it, recursively (`v_value(value(x))`), so that nested numbers (`ForwardDiff.hessian`) lose every level: stripping only the outer one leaks the inner derivatives into the payoff argument and gives wrong second derivatives. The base cases are `AbstractFloat` and `Integer` (an `EuropeanOption(1, 100)` keeps an `Int` maturity).
- **`v_mod(x)`** has a generic `exp`-based fallback. Backends may specialize it for speed or exactness (Dual, Hyper, only over `AbstractFloat` values; nested numbers use the fallback).
- **`lrm_interface_vec!(mu, sigma, args, mcBaseData, vb_mc)`** does the inner integration. It dispatches on the **MC config type × integration method**:
- **GPU** (`FinancialMonteCarlo.CudaMode`, with `using CUDA`): no extension is needed, `AbstractCudaMode`, `CudaMonteCarloConfig` and `get_matrix_type`/`get_array_type` exist without CUDA, and every transfer is a copy of a whole array (`to_host`, `on_device`), so no GPU array is indexed by scalars.
  - The metrics call `path_mean(lrm_paths(mu, sigma, args, mcBaseData, vb_mc))`. `lrm_paths` broadcasts `lrm_interface_vec!` per path; for a `CudaConfig` it evaluates the integrand on a matrix of paths × nodes in one broadcast and reduces along the nodes (`device_samples`, `device_option_paths`), since a kernel can run neither the loops of `lrm_interface_vec!` nor take the configuration (CPU rng). A smile is a matrix with one column per option (`path_mean`, `subordinator_mean`). `VibratoMonteCarloStandard` draws its inner samples with CURAND.
  - On the device: the paths of FinancialMonteCarlo, `simulate_with_vol` (Heston, same code as the CPU with `normals!`), the integration. On the host, then uploaded: `G_T` (Distributions has no GPU sampler) and the Kou jumps (path by path with the CPU rng). On the host entirely: the FFT density (a spline per path), `Vibrato` (regressions, backward induction), `PureLRM` (one point).
  - AD numbers must be isbits on the device: hyper duals and ForwardDiff (nested too) are tested, DualNumbers is isbits but untested, TaylorSeries and ReverseDiff cannot run there.
  - With `VibratoMonteCarloStandard` the inner samples are drawn again for every outer path from `mcBaseData.parallelMode.rng` (the `Z` from `init_lrm_vec` is only a buffer). Sharing them across paths turns their error into a bias that does not average out (deltas were 8-15% off).
  - `SerialMonteCarloConfig` × `VibratoMonteCarloStandard`: plain MC average.
  - `SerialAntitheticMonteCarloConfig` × `Standard`: antithetic average.
  - `SerialSobolMonteCarloConfig` × `Standard`: `cld(Nsim_lrm, 2)` points of a 1-D Sobol sequence, randomly shifted modulo 1 per path (Cranley-Patterson), through `norminvcdf`, each paired with its reflection.
  - any config × `VibratoMonteCarloAnalytic`: quadrature on a uniform grid `[x_min, x_max]` of `Npoints` (default 200 on ±7), always in standard deviations of the log increment around its mean. On Gaussian paths that is the standard normal variable. In `lrm_interface_aug!` (FFT density, non-Ito processes) the integral is a trapezoid with `PatternVector` weights, and the mean and std come from the FFT density (`DensitySpline`, built by `spline_density`); the grid must stay inside the spline domain (±15 in log-return for `PureLRM`, ±18 for `Vibrato`).
  - `VibratoMonteCarloGaussLegendre` (exported, `Npoints = 48` per piece): Gauss–Legendre on `[-8.5, 8.5 + σ]` standard deviations (`GAUSS_LEGENDRE_Z`; the `+ σ` covers the `e^{σz}` growth of a payoff linear in `S`), split at the payoff's singular points (`singular_points(option)`: `log K` for European, binary and American options, every strike for smiles, none for other payoffs). Each piece is smooth, so it is exact to about 1e-14 (1e-10 for call volgas at `σ√τ` ≈ 0.45), where Gauss–Hermite converges slowly across the kink or jump at the strike (13% of a binary's value per path with 200 nodes). It covers every `lrm_interface_vec!` caller, not the FFT density path (`lrm_interface_aug!`).
  - `VibratoMonteCarloGaussHermite`: Gauss–Hermite. It is not exported.
  - Both Gauss rules take their nodes from FastGaussQuadrature (`gausshermite`, `gausslegendre`), a regular dependency.

### Integration-method types ([src/utils.jl](src/utils.jl))
- `AbstractVibrato` is the supertype.
- Each type implements `init_lrm_vec(vb_mc, mcBaseData)`, which returns the inner sample points or grid.
- `Base.broadcastable` returns `Ref`, so these types are never broadcast over.
- `LRMInterfaceArguments(Z, option)` and `IntegrandArguments(mu, sigma, option)` bundle the vector and option arguments so that `@.` broadcasts only over `mu_jump`/`sigma_jump`. Keep this pattern when adding methods.

### Public metrics (dispatch on process type)
The package plugs into FinancialMonteCarlo's interface, `pricer(mcProcess, rfCurve, method, payoff)`, as FinancialFFT's `LewisMethod` does: each estimator is a method type in [src/methods.jl](src/methods.jl), subtype of `FinancialMonteCarlo.AbstractMethod` (through `AbstractVibratoMethod`), holding the `MonteCarloConfiguration` (`mc`) and, where there is one, the inner integration (`vb_mc`). The free functions `vibrato`, `vibrato_saltando`, `vibrato_am`, `pure_lrm` were removed after 0.2.0. European methods take `EuropeanLike` (an `EuropeanPayoff` or an `AbstractSmile`), so that `Vibrato` on an `AmericanPayoff` dispatches to the American method without ambiguities. Arrays of payoffs, `Position` and `Portfolio` are not supported (FinancialMonteCarlo types those `pricer` methods on `MonteCarloConfiguration`).
- [vibrato.jl](src/metrics/vibrato.jl), `pricer(process, rfCurve, Vibrato(mc, vb_mc), option)`:
  - `BlackScholesProcess` / `ItoProcess`: Gaussian last step.
  - `HestonProcess`: uses `simulate_with_last_vol` to get the last-step variance. Serial, Antithetic and Sobol configs (`simulate_with_vol` has one method each: FinancialMonteCarlo's `SerialMonteCarloConfig` matches any method, so a missing one silently falls back to the rng), `ZeroRate` only.
  - any other process: FFT-derived spline density of the last step (`spline_density` + `lrm_interface_aug!` in `lrm_aug.jl`). This path covers Lévy processes.
  - `VibratoMonteCarlo.Conditional(mc)` (not exported) is an analytic-pricer baseline (BS only).
- [vibrato_saltando.jl](src/metrics/vibrato_saltando.jl), `pricer(process, rfCurve, VibratoSaltando(mc, vb_mc), option)`:
  - `FiniteActivityProcess` (Kou, Merton): conditions on the **last jump time** (`simulate_log_with_jump_times`) and integrates the Gaussian diffusion from there to `T`. Jumps are applied at the next grid point, so the simulator returns the grid point *before* the step containing the last jump (plus that step's jumps). This keeps the conditional step exact and at least one `dt` long; rounding up instead gives `T - T_jump = 0` and NaN derivatives. Differentiating w.r.t. `λ` or `T` is unsupported and only emits a warning.
  - `InfiniteActivityProcess` (VG, NIG): conditions on the whole subordinator path, i.e. on `G_T ~ FinancialMonteCarlo.compute_subordinator_law(process, κ, T)`. Given `G_T`, `log S_T` is exactly Gaussian with mean `log S0 + (r - d - ψ) T + θ G_T` and variance `σ² G_T`, so `Nstep` does not matter. Conditioning on the last step only is exact too, but for VG that increment (Gamma, shape `dt / κ`) is often so small that `σ √ΔG` is lost to rounding against `log S`. `G_T` comes from `simulate_subordinator` in `src/simulation/subordinator.jl`; the path simulators of the subordinated Brownian motion were removed, FinancialMonteCarlo has them.
    - `κ` and `T` are supported, unlike `λ`: they change the law of `G_T`, and a draw is not a smooth function of them (rejection samplers). `simulate_subordinator` draws `G_T` from their values and returns the weight `v_mod(subordinator_log_density(...))`, the LRM step applied to the outer law. Without derivatives on them (`carries_derivatives`) the weight is `true`, so no density is computed.
    - `subordinator_mean` takes out the sample mean as a baseline, which cuts the noise of the `κ` and `T` derivatives 3-19×. Its input must be discounted already: a discount outside the baseline multiplies its `T` derivative by the score noise.
    - `log_gamma` is elementary because DualNumbers, HyperDualNumbers and TaylorSeries have no `loggamma`.
  - `VibratoMonteCarlo.ConditionalSaltando(mc)` (not exported): same paths as `VibratoSaltando`, but prices what the conditioning leaves in closed form (`analytic_pricer`, European and binary only) and differentiates pathwise, with the same weight on the law of `G_T` for infinite activity. After the last jump for `FiniteActivityProcess`; for `InfiniteActivityProcess` Black-Scholes with volatility `σ √(G_T / T)` and the drift as a dividend yield. It is the exact inner-integration limit of `VibratoSaltando`. For infinite activity the price is then a 1-D integral over `G_T` (Madan-Carr-Chang for VG), so it is a test oracle more than a pricer: Lewis is the better European pricer.
- [pure_lrm.jl](src/metrics/pure_lrm.jl), `pricer(process, rfCurve, PureLRM(mc, vb_mc), option)`: no path simulation. Pure LRM over the full-maturity density: Gaussian for `ItoProcess`, FFT spline otherwise.
- [vibrato_am.jl](src/metrics/vibrato_am.jl), `pricer(process, rfCurve, Vibrato(mc, vb_mc), am_opt::AmericanPayoff)`: exercise at the `Nstep` grid dates, policy from Longstaff–Schwartz (regression on `1, S, S²` of the in-the-money paths, `S` scaled by `S0`, `v_value`s only). Two methods:
  - `BlackScholesProcess`, vibrato over every exercise date. Pass 1 (independent paths, values only, `exercise_regions`) turns the regressions into exercise regions in log price (`positive_intervals`), fixed in price space. Pass 2 simulates `2 Nstep` half steps (same MC configuration, shared rng) so that every exercise date `t_k` is the midpoint of `[s_k, s_{k+1}]`, `s_k = (k - 1/2) dt`, where `log S_{t_k}` is a Brownian bridge (mean of the ends, std `σ √dt / 2`). Backwards along each path, `bridge_vibrato` integrates "exercise if in the region, else continue with the path value" against the bridge density (`lrm_segment`: Gauss–Legendre on each region piece, points from values, derivatives in `v_mod`); the payoff at `T` gets the usual vibrato over the last half step. Since the optimal boundary does not depend on `S0`, a policy fixed in price space gives the `S0` derivatives of every order (the regression error enters them only at second order).
  - Lévy and Heston, frozen decisions. The decision to exercise at `t_k` is taken at `t_{k-1}` and the cash flow is the vibrato expectation `E[payoff(S_{t_k}) | S_{t_{k-1}}]` over one step (`frozen_laws`: FFT `DensitySpline` for Lévy, which needs `VibratoMonteCarloAnalytic`; for Heston a `GaussianStep` with the variance from `simulate_with_vol`). One AD cash flow per path at its exercise date. Other `ItoProcess`es throw.

### Options
- `EuropeanOption`, `BinaryEuropeanOption`, or a smile of them, `EuropeanOptionSmile` or `BinaryEuropeanOptionSmile` (FinancialMonteCarlo, not exported), which has **one** maturity and returns a vector of prices. Arrays of options are not supported.
- Never reimplement payoffs: `payout_f` on a smile calls `FinancialMonteCarlo.payout` through `smile_to_options`, as FinancialMonteCarlo's own `pricer` does. That allocates a vector of options at every node; the proper `payout(S, smile)` is a TODO of FinancialMonteCarlo, tracked in its CLAUDE.md.
- `Vibrato` on American options takes one FinancialMonteCarlo `AmericanPayoff` (`AmericanOption`, `BinaryAmericanOption`); the integrands call `FinancialMonteCarlo.payout`, so any payoff that defines it works.
- `analytic_pricer` exists only for European and Binary European options.

### Simulation helpers ([src/simulation/](src/simulation/))
These are custom simulators that return extra state the metrics need: the last jump time, the variance path (`simulate_with_vol`; `simulate_with_last_vol` keeps its last column), or the subordinator `G_T` with its weight (`simulate_subordinator`). They allocate matrices with an element type derived from a "zero" built out of all parameters (e.g. `isDualZero`), so AD numbers propagate. Follow this pattern in new simulators. Use `eps_value(x)`, not `eps(typeof(x))`: TaylorSeries does not define `eps`.

### AD backend glue: package extensions ([ext/](ext/))
- Each AD backend is a Julia package extension (`VibratoMonteCarlo<Pkg>Ext`). The backends are declared under `[weakdeps]`/`[extensions]` in `Project.toml`, which requires Julia ≥ 1.9 (the package needs 1.10).
- To support a new number type:
  1. Add `ext/VibratoMonteCarlo<Pkg>Ext.jl` that does `import VibratoMonteCarlo: v_value, v_mod` and adds methods.
  2. Register the package in `[weakdeps]`, `[extensions]` and `[compat]`.
  3. Add it to `test/Project.toml`.
- Extensions cannot define new types in the parent module, nor replace a method the parent already defines (so an optional dependency cannot swap in a better implementation when loaded).
- Mooncake differs from the other backends: it sees `Float64` inputs, not a number type, so `VibratoMonteCarloMooncakeExt` adds rules instead of `v_value`/`v_mod` methods. `v_value` is a zero-derivative primitive and `v_mod` is overlaid with its `exp` form, because `v_mod(::AbstractFloat) = 1` would drop every score term (split Gauss-Legendre would still be exact through its moving cut, the other rules would lose the jump of a binary). `carries_derivatives` is overlaid to `true` on floats too, otherwise `κ` and `T` would lose the weight on the law of `G_T`. `exercise_regions` is zero-derivative (foreigncall) (`set_seed!` has its rule in FinancialMonteCarlo's own Mooncake extension). The rule for `fft` (FFTW plans use `try/finally`) comes from FinancialFFT's Mooncake extension, which has no `bfft` rule: forward over reverse through the FFT density would need one. Use `Mooncake.Config(enable_nfwd = false)` in forward mode: the default path pushes `NDual` numbers through the code and fails. Known gaps: forward mode on arrays of options, since removed, did not finish in 2 h (reverse worked), smiles were not measured, reverse and Hessians of `Vibrato` on Lévy and Hessians of the FFT `Vibrato` fail on `copy(::TwicePrecision)` in Interpolations (forward works), the Hessian with `VibratoMonteCarloStandard` fails. The first call compiles for 2-30 min on each case, so ForwardDiff stays the better choice for few parameters.
- ForwardDiff is pulled in transitively (FinancialFFT 0.3 → DataInterpolations), so its extension is usually active even if the user never loads it. FinancialFFT `master` drops DataInterpolations, so after the next FinancialFFT release this no longer holds.
- `src/deps/zygote_dependencies.jl` is untracked and not wired in (marked `#wrong`).

## Untracked and experimental files — do not treat as package code

`src/VibratoMonteCarlo.jl` does **not** `include` the following files in `src/metrics/`. They are scratch scripts (many are untracked, some are duplicate "copy"/"Copia" files):
- `vibrato_path.jl`
- `vibrato_zygote*.jl`, `simulate_zygote*.jl`: Zygote reverse-mode experiments with hand-written `rrule`s.

Don't wire them into the module, and don't delete them, unless asked.

## Conventions

- Formatting: JuliaFormatter with `margin = 1200`, so long signatures stay on one line. Also `whitespace_typedefs = true` and `trailing_comma = false`.
- Every public metric calls `FinancialMonteCarlo.set_seed!(mcBaseData)` first, so results are reproducible. Keep this in new metrics.
- Reuse FinancialMonteCarlo, FinancialFFT and FinancialToolbox before writing a helper: payoffs (`payout`), laws (`compute_subordinator_law`), closed forms (`blsprice`, `blsbin`). Write one only when they have nothing generic in the number type (`log_gamma`, `subordinator_log_density`: `loggamma` and `Distributions.logpdf` are not).
- Drifts come from FinancialFFT: `drift_correction(process) = FinancialFFT.characteristic_exponent_i(1, process)`, i.e. `log E[e^{X_1}]` of the driftless log price, and the risk neutral drift of `log S` is `r - d - drift_correction`. Never write `σ² / 2` or `λ * compute_drift` by hand. A process without the two-argument `characteristic_exponent_i` (Heston, the log-normal mixtures) throws in the generic `ItoProcess` methods.
- Write `a * b + c` under `@muladd` (or `muladd`) in new code.
- The code uses `@muladd` (MuladdMacro), `@views`, and broadcast-heavy `@.` style. Constants come from IrrationalConstants (`log2π`, `sqrt2`, `sqrtπ`).
- Prices are discounted with `exp(-r * T)` at the end of each metric.

## CI and landing page

- **GitHub Actions** is the CI:
  - `CI.yml` runs the test suite on every push (Julia 1.10 and latest), with the private registry added through `julia-buildpkg`'s `localregistry` input, and uploads the coverage to Codecov (needs the `CODECOV_TOKEN` secret).
  - `CompatHelper.yml` runs nightly and opens the compat-bump PRs that make up most recent commits.
  - Every action runs on Node.js 24 (`checkout@v7`, `setup-julia@v3`, `cache@v3`, `codecov-action@v7`, the Pages actions; the `julia-buildpkg`, `julia-runtest` and `julia-processcoverage` v1 are composite). CompatHelper installs Julia with `setup-julia` instead of relying on the runner image. Keep new actions on Node 24 versions.
  - `DeployPage.yml` builds the Franklin/PkgPage site in `page/` (math write-up in `page/index.md`) to `gh-pages`. It copies `etc/logo.png` to `page/_assets/logo.png` first, for the navbar and the favicon: that copy is gitignored, the logo lives in `etc/` only.
- GitLab is no longer used: the `origin2` remote is legacy. Push only to `origin` (GitHub).
- The logo is drawn by `etc/logo.jl` (Luxor and the Georgia font, not in any Project.toml).

The landing page (`page/index.md`):
- Derived formulas are never typed by hand: `#hideall` Julia blocks generate them with Symbolics at build time from `VibratoMonteCarlo.log_density`, using the helpers of the first blocks (`latex_weight`, `generate_faa_di_bruno_exp`), and `\textoutput{}` shows them. Simplify with `δ = √h` as a plain symbol, since Symbolics does not cancel `sqrt(T)`. Definitions, model inputs and general identities stay in the text. Check new formulas numerically before publishing them.
- Blocks run in page order, so a helper must come before its first use, and a Julia string needs `\\` for every LaTeX backslash.
- Subsections fold with `\begin{expand}{title}` (`config.md`, `<details>`, style in `_css/custom.css`). Keep the summary table and the section intros open.
- Citations are `\cite{key}` with `\biblabel{key}{Name (Year)}` in the References section.
- To preview, copy `page/` to a scratch folder, `Pkg.develop` this package there, copy the logo, and run `Franklin.optimize(prerender = false, minify = false)` (about 3 minutes on Julia 1.10). `Franklin.serve` hung here. The output in `__site` uses the `VibratoMonteCarlo.jl` prepath, so serve its parent with a static server.
