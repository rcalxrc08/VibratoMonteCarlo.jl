# Tests

_What each test file covers, seeds and tolerances, durations._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

`test/runtests.jl` includes an **explicit list** of files, so new test files must be added there. Each file guards on
`TEST_UTILS_LOADED`, so you can also `include` it on its own. The suite takes about 12 min on Julia 1.11; the GPU tests
about 4 min on an MX150.

## Files

- `utils.jl`: loads every AD backend, defines model constructors (`kou`, `merton`, `heston`, `variance_gamma`) and
  `test_value_and_derivatives` (value, 1st and 2nd derivative vs a reference, for `Hyper` and `Taylor1`).
- `context.jl`: the context layers follow the code (every file a contract lists exists, every included source file and
  every reference is reached from a contract or CLAUDE.md, every contract is routed from the root CONTEXT.md). Seconds.
- `backends.jl`: all five AD extensions give the same BS price, delta and gamma, including nested ForwardDiff
  (`ForwardDiff.hessian`). No method ambiguities in `pricer` (`Test.detect_ambiguities`).
- `vibrato_models.jl`: Kou/BS/Heston/VG × {standard+analytic, standard+MC, antithetic+MC, Sobol+MC}. Compares against
  Lewis, `blsprice`, or the MC pricer (Heston). `PureLRM` on Heston through the FFT density against Lewis (calls at three
  strikes; the gradient in `S0` and the vol of vol on ±10 standard deviations, 1e-5).
- `forwarddiff.jl`: multi-partial ForwardDiff duals (Merton, BS).
- `taylor.jl`: `Taylor1` inputs with Gauss–Hermite integration, against references; on Heston and the FFT density
  against hyper duals on the same paths; `eps_value`; `TaylorN` on Heston against ForwardDiff; `TaylorN` and
  `ForwardDiff.hessian` through the FFT density (Kou) against hyper duals.
- `regressions.jl`: one case per fixed bug (Heston antithetic correlation sign, Sobol ignored on Heston and on `G_T`,
  antithetic ignored on `G_T` (#29), Sobol and antithetic ignored on the Kou and Merton jumps (#30), `PureLRM` default
  Analytic grid, Kou last-jump-time alignment, shared inner samples, ±5σ default Analytic bounds). The antithetic Kou
  test (std of the delta, antithetic against plain) runs 256 seeds with a threshold of 1.2: the true ratio is about 1.45
  (1.13: 1.51 on 384 seeds; 1.10, whose MersenneTwister seeding differs: 1.28 and 1.37 on two blocks of 128, 1.44 before
  FinancialMonteCarlo#54), ±0.13 on 128 seeds; on seeds 1:32 it went from 1.78 to 1.37 with FinancialMonteCarlo 0.9.
- `levy.jl`: VG and NIG value, delta and gamma against Lewis, for `VibratoSaltando` and the FFT `Vibrato`. Derivatives
  over `κ` and `T`: `log_gamma`, `subordinator_log_density` against `logpdf`, `ConditionalSaltando` against Lewis, all
  AD backends agree, a smile against its options alone, and no weight without derivatives.
- `american.jl`: `Vibrato` on American options. BS put value, delta and gamma against a Crank–Nicolson PDE of the
  Bermudan with the same dates (`bermudan_pde`, defined in the file); calls without dividends against the European price
  (BS and Kou with gamma, Heston value and delta against the MC pricer); all AD backends agree.
- `conditional.jl`: `ConditionalSaltando` (Kou, Merton, VG, NIG; European and binary) against `VibratoSaltando` on the
  same paths (Gauss–Hermite within its quadrature error, Gauss–Legendre to 1e-12, on VG and NIG also over `κ` and `T`),
  `Conditional` against `Vibrato` on BS, and the warning when differentiating over `λ`.
- `gauss_legendre.jl`: `VibratoMonteCarloGaussLegendre` through `PureLRM` on BS against `blsprice`/`blsbin` to 1e-9
  (calls, puts, binaries, short and long maturities, strikes outside the bounds), with ForwardDiff Hessians, ReverseDiff,
  DualNumbers, `Taylor1` and a smile.
- `smiles.jl`: European and binary smiles through `Vibrato`, `VibratoSaltando` and `PureLRM`: every entry against the
  option priced alone on the same seed, a binary smile against `blsbin`.
- `mooncake.jl`: reverse gradients of a BS binary and of `VibratoSaltando` on VG (with `κ` and `T`) against ForwardDiff,
  forward gradient of the FFT `Vibrato` on VG against hyper duals (ForwardDiff against them too). Compiling the rules
  takes minutes, so it has few paths and few cases.
- `check/accuracy.jl`: tight `1e-3` checks against Lewis. **Opt-in** (`VIBRATO_TEST_ACCURACY=true`). It currently fails
  (`references/method_limits.md`).
- `examples.jl`: README.md is up to date with the examples (`etc/readme.jl` copies `examples/<name>.jl` between
  `<!-- BEGIN examples/<name>.jl -->` and `<!-- END ... -->`: GitHub cannot include a file, so after changing an example
  run `julia etc/readme.jl`), and `examples/hello_world.jl` runs and matches `blsprice`.
- `quality.jl`: Aqua (`Aqua.test_all`), `@inferred` on every pricer family, the execution mode of the methods
  (`execution_mode`, `with_mode`, batched modes refused), the input validation (`ArgumentError` from the
  inner-integration constructors, `DomainError` for `T <= 0` in the Heston simulations), and Float32 prices from Float32
  models.
- `cuda/runtests.jl`: GPU tests, **opt-in** and in their own project (`test/cuda/Project.toml`, with `[sources]` for the
  package, Julia ≥ 1.11), so CI neither downloads nor compiles CUDA. With `CUDA.allowscalar(false)`: the device inner
  integration against the host one on the same `mu` (≤ 4e-11), the GPU Heston paths against FinancialMonteCarlo's GPU
  Heston (1.4e-8, the variance floors differ), `G_T` of plain MC drawn on the device (its distribution function at the
  deciles within 4 standard errors, reproducible with `set_seed!`, the weight of `κ` and `T` against the host on the same
  draws), antithetic VG/NIG on the same `G_T` as the host (1e-15, halves paired as `path_sides`), the last jumps on the
  device against the old per path algorithm (the jumps of the same `JumpList` added to a host copy of each row in time
  order, 1e-15), the FFT density on the device (the FFTs within 2.1e-15, the integration of the same paths on the same
  law within 1.2e-13 on values and 1.1e-12 on gammas, `PureLRM` against the host to 1e-9, a nested ForwardDiff Hessian
  of the FFT `Vibrato` against hyper duals, Kou and NIG against Lewis), the CURAND estimators against references
  (tolerances twice the worst of 6 seeds, VG and NIG with the device `G_T` among them), nested ForwardDiff on the device
  against hyper duals, American options (`bridge_step` and `step_payoffs` against the host on the same paths and laws,
  the BS and Heston put against the host pricer, the Kou call without dividends against the European Lewis price) and
  `PureLRM` against the host (Black-Scholes 1e-15 with a quadrature, Heston through the FFT density 4.2e-10, with a
  nested ForwardDiff Hessian).
- `cuda/timings.jl`: not a test, the GPU timings of the parts moved by #51 (`references/gpu.md`).
- `runner_z.jl` is not included: it needs `vibrato_z` from the untracked Zygote experiments.

## Seeds and tolerances

The suite uses fixed seeds (`set_seed!`), so results are deterministic for a given Julia version. The seeded
MersenneTwister stream differs between Julia 1.10 and 1.11+, so calibrate a tolerance on the spread over several seeds
(`MonteCarloConfiguration(Nsim, Nstep; seed = seed)`), not on one run. If a numerics change trips one, check the error
size before loosening it.

Known package bugs are pinned with `@test_broken`, which turns into an error once a bug is fixed. That is the signal to
change it to `@test`. There are none at the moment.
