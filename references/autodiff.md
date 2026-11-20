# Automatic differentiation: number types, extensions, Mooncake

_Read before touching `v_value`, `v_mod`, an extension, or the element types of a simulation._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

## Extensions (`ext/`)

- Each AD backend is a Julia package extension (`VibratoMonteCarlo<Pkg>Ext`): DualNumbers, HyperDualNumbers,
  ForwardDiff, ReverseDiff, TaylorSeries, Mooncake. The backends are declared under `[weakdeps]`/`[extensions]` in
  `Project.toml`, which requires Julia ≥ 1.9 (the package needs 1.10).
- To support a new number type:
  1. Add `ext/VibratoMonteCarlo<Pkg>Ext.jl` that does `import VibratoMonteCarlo: v_value, v_mod` and adds methods
     (`v_value` recursive, `references/inner_integration.md`).
  2. Register the package in `[weakdeps]`, `[extensions]` and `[compat]`.
  3. Add it to `test/Project.toml`.
- Extensions cannot define new types in the parent module, nor replace a method the parent already defines (so an
  optional dependency cannot swap in a better implementation when loaded).
- ForwardDiff is not a dependency (FinancialFFT 0.4 dropped DataInterpolations, which pulled it in): its extension is
  active only when the user loads it.

## Mooncake

Mooncake differs from the other backends: it sees `Float64` inputs, not a number type, so
`VibratoMonteCarloMooncakeExt` adds rules instead of `v_value`/`v_mod` methods.

- `v_value` is a zero-derivative primitive and `v_mod` is overlaid with its `exp` form, because
  `v_mod(::AbstractFloat) = 1` would drop every score term (split Gauss-Legendre would still be exact through its moving
  cut, the other rules would lose the jump of a binary).
- `carries_derivatives` is overlaid to `true` on floats too, otherwise `κ` and `T` would lose the weight on the law of
  `G_T`.
- `exercise_regions` is zero-derivative (foreigncall) (`set_seed!` has its rule in FinancialMonteCarlo's own Mooncake
  extension). The rules for the FFT of the density come from FinancialFFT's Mooncake extension.
- Use `Mooncake.Config(enable_nfwd = false)` in forward mode: the default path pushes `NDual` numbers through the code
  and fails.
- Known gaps: forward mode on arrays of options, since removed, did not finish in 2 h (reverse worked), smiles were not
  measured, reverse and Hessians of `Vibrato` on American options on Lévy and Hessians of the FFT `Vibrato` failed on
  `copy(::TwicePrecision)` in Interpolations (forward worked), measured before the switch to FinancialFFT 0.4 and not
  measured again, the Hessian with `VibratoMonteCarloStandard` fails.
- The first call compiles for 2-30 min on each case, so ForwardDiff stays the better choice for few parameters.

## Number types (#45)

- `init_lrm_vec(vb_mc, mcBaseData, paths)` builds the inner samples, nodes and grids in `rng_type(paths)`
  (FinancialMonteCarlo's `get_rng_type` of the element type: Float32 for Float32 paths, Float64 for AD numbers), as do
  `node_row` and the CURAND samples on the GPU, and the Gauss-Legendre bounds take the type of `sigma`; Float64 models
  draw exactly the same numbers as before.
- A Float32 model gives a Float32 price (`test/quality.jl`), except `VibratoSaltando` on VG/NIG (`draw_subordinator`
  fills a Float64 buffer) and American options.
- Payoffs and smiles broadcast as scalars through FinancialMonteCarlo's `Base.broadcastable`: no `Ref(option)`.
- Simulators derive their element type from a zero built out of all the parameters (`references/simulation.md`).

## Known gaps of the number types

- The FFT density path (`Vibrato` and `PureLRM` on non-Ito processes, European and American) takes `TaylorN` and nested
  ForwardDiff duals (`ForwardDiff.hessian`) (`[compat]` requires FinancialFFT 0.6, whose `FFTDensityInverter` is
  generic in the number type since 0.4). Hessians match hyper duals to ≈ 1e-10 relative on VG, NIG, Kou and Merton.
- `Vibrato` through the FFT density (European and American) does not take derivatives over `κ` on VG and NIG, for any
  number type: the path simulation builds `Gamma` with AD parameters. `PureLRM` does, and so does `VibratoSaltando`.
- `VibratoSaltando` on Kou and Merton: derivatives over `λ` or `T` are unsupported (a warning).
- On the GPU only isbits numbers: `references/gpu.md`.

## Zygote

`src/deps/zygote_dependencies.jl` is untracked and not wired in (marked `#wrong`); `test/runner_z.jl` needs `vibrato_z`
from the untracked Zygote experiments and is not included in the suite.
