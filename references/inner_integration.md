# The vibrato step and the inner integrations

_Stage 01 and 03 reference: the core idea, `v_value`/`v_mod`, the rules of the inner integration._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

## Core idea

1. Simulate the path with ordinary Monte Carlo up to `T - dt` (or up to the last jump time).
2. Price the last step by integrating the payoff against the known transition density, using the **likelihood ratio
   method (LRM)**. Sample points are built from `v_value(mu) + v_value(sigma) * z`, so the payoff argument carries no
   derivative. All parameter dependence flows through the density weight.
3. The density weight is `v_mod(log_density)`, where `v_mod(x) = exp(x - v_value(x))`. It has value 1, and its
   derivatives are the score terms at every order (an automatic Faà di Bruno).
4. The integrand is `(payoff - baseline) * v_mod(log_density) + baseline`, with `baseline` the payoff at the centre of
   the density (no derivative). The value is unchanged; the score terms have zero mean, so the derivatives lose the
   baseline times their integration error. This is the vibrato control variate: it cuts the MC-inner gamma noise 10-80×.

All of this sits in `src/03_inner_integration/lrm.jl`:

- **`v_value(x)`** strips the derivative part off an AD number. Each AD backend extension **must** implement it,
  recursively (`v_value(value(x))`), so that nested numbers (`ForwardDiff.hessian`) lose every level: stripping only the
  outer one leaks the inner derivatives into the payoff argument and gives wrong second derivatives. The base cases are
  `AbstractFloat` and `Integer` (an `EuropeanOption(1, 100)` keeps an `Int` maturity).
- **`v_mod(x)`** has a generic `exp`-based fallback. Backends may specialize it for speed or exactness (Dual, Hyper, only
  over `AbstractFloat` values; nested numbers use the fallback).
- **`lrm_interface_vec!(mu, sigma, args, mcBaseData, vb_mc)`** does the inner integration. It dispatches on the **MC
  config type × integration method** (below).
- `eps_value(x)`, not `eps(typeof(x))`: TaylorSeries does not define `eps`. `carries_derivatives(x)` is false on floats
  and integers.

## Integration-method types (`src/01_methods/utils.jl`)

- `AbstractVibrato` is the supertype.
- Each type implements `init_lrm_vec(vb_mc, mcBaseData, paths)`, which returns the inner sample points or grid.
- `Base.broadcastable` returns `Ref`, so these types are never broadcast over.
- `LRMInterfaceArguments(Z, option)` (`src/04_european/vibrato.jl`) and `IntegrandArguments(mu, sigma, option)`
  (`lrm.jl`) bundle the vector and option arguments so that `@.` broadcasts only over `mu_jump`/`sigma_jump`. Keep this
  pattern when adding methods.

## The rules, config × inner integration

- `MonteCarloConfiguration` × `VibratoMonteCarloStandard`: plain MC average.
- `MonteCarloConfiguration{<:AntitheticMC}` × `Standard`: antithetic average.
- `SobolConfig` × `Standard`: `cld(Nsim_lrm, 2)` points of a 1-D Sobol sequence, randomly shifted modulo 1 per path
  (Cranley-Patterson), through `norminvcdf`, each paired with its reflection.
- With `VibratoMonteCarloStandard` the inner samples are drawn again for every outer path from
  `FinancialMonteCarlo.host_rng(mcBaseData)` (the `Z` from `init_lrm_vec` is only a buffer). Sharing them across paths
  turns their error into a bias that does not average out (deltas were 8-15% off).
- any config × `VibratoMonteCarloAnalytic`: the trapezoid rule on a uniform grid `[x_min, x_max]` of `Npoints` (default
  200 on ±7), always in standard deviations of the log increment around its mean, with the weights 1, 2, ..., 2, 1 of
  `trapezoid_weights` (a `PatternVector`) times `step / 2`, the step read from the range (it was the rectangle rule, with
  the step from two reductions, until #69). `init_lrm_vec` returns the nodes and these weights times the normal density,
  computed once per price: an `exp` per node and path made the rule twice as slow (#74; 16 ms for a binary on 10 000
  hyper-dual paths with 20 nodes). The FFT density path takes the range itself (`analytic_grid`). On Gaussian paths that is the standard normal variable. In
  `lrm_interface_aug!` (`lrm_aug.jl`: FFT density, non-Ito processes) the integral is the same trapezoid over the nodes
  of one path, evaluated with one call of the vectorized
  `FinancialFFT.density_values` and a broadcast of the integrand; the mean and std come from the FFT density
  (`FFTStepLaw`, built by `fft_step_law`). Outside the Chernoff range of the density its value is zero, floored at `eps`
  (`floor_density`).
- `VibratoMonteCarloGaussLegendre` (exported, `Npoints = 48` per piece): Gauss–Legendre on `[-8.5, 8.5 + σ]` standard
  deviations (`GAUSS_LEGENDRE_Z`; the `+ σ` covers the `e^{σz}` growth of a payoff linear in `S`), split at the payoff's
  singular points (`singular_points(option)`: `log K` for European, binary and American options, every strike for
  smiles, none for other payoffs). Each piece is smooth, so it is exact to about 1e-14 (1e-10 for call volgas at `σ√τ` ≈
  0.45), where Gauss–Hermite converges slowly across the kink or jump at the strike (13% of a binary's value per path
  with 200 nodes). It covers every `lrm_interface_vec!` caller, not the FFT density path (`lrm_interface_aug!`).
- `VibratoMonteCarloGaussHermite`: Gauss–Hermite. It is not exported.
- Both Gauss rules take their nodes from FastGaussQuadrature (`gausshermite`, `gausslegendre`), a regular dependency.

## Smiles and inference

A smile is integrated strike by strike on every path (`map_strikes`): each strike is a `SmileStrike{O}(K, isCall)`,
whose payoff is FinancialMonteCarlo's `payout_strike`, the formula a smile shares with its single options, so no vector
of options is built per path (`smile_to_options` is not used). Fixed-node rules (Gauss-Legendre, Analytic,
Gauss-Hermite) integrate each strike on the same nodes (`vibrato.jl`, #42: the integrand of the whole smile built its
options and two vectors at every node, 182 MB for 2000 paths); `VibratoMonteCarloStandard` draws its samples once per
path and integrates each strike on them (`inner_samples!`, `inner_mean`), the FFT density path computes the density of
the path once and sums each strike on it (`trapezoid_sum`). The sums take the integrand at each node without an array of
its values. Allocations per path, 1000 paths, Float64 and hyper duals alike (#74): a single option 2 → 0, a smile of two
with `Standard` 227-447 → 4, with the fixed-node rules 6-10 → 4, on the FFT density (VG) 293 → 46; prices unchanged to
2e-15. Gauss-Legendre sums its pieces with `map` over
`segment_starts`/`segment_ends` (a tuple for an option, a vector for a smile), not a filtered generator, which is not
inferred: zero-length pieces contribute zero. Every pricer is inferred (`test/quality.jl`, with Aqua), on Julia 1.10
too: there a generator whose elements build a closure is not inferred (see `references/american.md` for the bridge).

## On the device

The metrics call `path_mean(lrm_paths(mu, sigma, args, mcBaseData, vb_mc))`. `lrm_paths` broadcasts
`lrm_interface_vec!` per path; for a `CudaConfig` it evaluates the integrand on a matrix of paths × nodes in one
broadcast and reduces along the nodes (`device_samples`, `device_option_paths`). Details in `references/gpu.md`.
