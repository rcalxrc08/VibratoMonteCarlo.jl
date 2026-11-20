# Dependencies and compat history

_What each of the author's packages provides, the registry, and what each FinancialMonteCarlo release changed here._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

## The author's packages

- **FinancialMonteCarlo.jl**: processes (`BlackScholesProcess`, `KouProcess`, `MertonProcess`, `HestonProcess`,
  `VarianceGammaProcess`, `NormalInverseGaussianProcess`), payoffs, `MonteCarloConfiguration`, `ZeroRate`, `simulate`,
  `pricer`, `set_seed!`, the jumps (`JumpList`, `add_jumps!`), the subordinator laws and draws, `draw_normals!`, modes.
- **FinancialFFT.jl**: the transition densities of the non-Ito processes (`density` with `FFTDensityInverter`, evaluated
  by `density_values`), `characteristic_exponent_i` for the drifts, and `LewisMethod` reference prices.
- **FinancialToolbox.jl** (0.8: `d` is a required argument of `blsprice`, `normcdf` uses `erfc`): `blsprice`/`blsbin`,
  used for closed-form checks and `analytic_pricer`.
- **PatternVectors.jl**: the trapezoid weights of the FFT density path (`lrm_aug.jl`).
- External: FastGaussQuadrature (`gausshermite`, `gausslegendre`), Distributions, Sobol, StatsFuns, MuladdMacro,
  IrrationalConstants. Interpolations and FFTW are no longer dependencies; ForwardDiff is a weak one
  (`references/autodiff.md`).

Some dependencies come from a private registry. Add it before `Pkg.instantiate()`:

```julia
import Pkg
Pkg.Registry.add(Pkg.RegistrySpec(url = "https://github.com/rcalxrc08/LocalRegistry.git"))
Pkg.Registry.add("General")
```

## FinancialMonteCarlo releases

- **0.9** (VibratoMonteCarlo 0.6.1; compat `"0.7, 0.8, 0.9"` in the package, tests, `test/cuda`, docs and `bench`): no
  change of the code. Its antithetic jumps no longer draw a uniform for the type of their table (FinancialMonteCarlo#54),
  which shifts the stream: the antithetic Kou test of `regressions.jl` was retuned (`references/testing.md`). Its curves
  are isbits and two curves only add or subtract (lazily); the drifts still assume scalar rates (#68).
- **0.8** (#63, VibratoMonteCarlo 0.5.1; compat `"0.7, 0.8"`): no change of the code was needed, the whole suite passes
  on it with FinancialFFT 0.6.1 (the only FFT release accepting 0.8). Its changes (invalid inputs throw
  `DomainError`/`DimensionMismatch`/`ArgumentError` instead of `AssertionError`, `ShiftedLogNormalMixture` removed,
  `LogNormalMixture` takes tuples) touch nothing used here. `vibrato_am.jl` rebuilds the configuration with the field
  constructor (`references/american.md`).
- **0.7 and FinancialFFT 0.6** (VibratoMonteCarlo 0.5.0): no change of the code was needed. The 0.7 breaking changes
  (`Position` as two vectors, tuple keys of `MarketDataSet`/`Portfolio`, `Curve` without Dictionaries/Interpolations,
  one-dimensional inputs as `AbstractVector`s) touch nothing used here; smiles still go through `smile_to_options`.
- **0.6** (#42): a `MonteCarloConfiguration{Sampling, Mode, State}` with fields `Nsim, Nstep, sampling, mode, rng,
  reseed`, built with keywords (`sampling`, `mode`, `seed`, `rng`, `reseed`); the host generator is
  `FinancialMonteCarlo.host_rng(mcBaseData)` (not `parallelMode.rng`), the configurations are dispatched as
  `MonteCarloConfiguration`, `MonteCarloConfiguration{<:AntitheticMC}` and `SobolConfig` (the old `Serial…Config`
  aliases are gone), and `get_matrix_type` takes a type. The batched modes are refused (`check_mode`).
