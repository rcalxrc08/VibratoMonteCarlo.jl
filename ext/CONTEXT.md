# Extensions

One extension per AD number type, plus the Mooncake rules.

## Inputs

| Load | Why |
|:--|:--|
| `VibratoMonteCarloDualNumbersExt.jl`, `VibratoMonteCarloHyperDualNumbersExt.jl`, `VibratoMonteCarloForwardDiffExt.jl`, `VibratoMonteCarloReverseDiffExt.jl`, `VibratoMonteCarloTaylorSeriesExt.jl` | `v_value` (recursive) and `v_mod` for each number type |
| `VibratoMonteCarloMooncakeExt.jl` | rules: `v_value` zero-derivative, `v_mod` as `exp`, `carries_derivatives` true, `exercise_regions` zero-derivative |
| `references/autodiff.md` | how to add a number type, the Mooncake gaps |

## Process

An extension only adds methods to `v_value`, `v_mod` (or Mooncake rules); types and functions live in `src/`. Register a
new one in `[weakdeps]`, `[extensions]`, `[compat]` and `test/Project.toml`.

## Outputs

Checked by `test/backends.jl` (every backend, nested ForwardDiff), `test/taylor.jl`, `test/mooncake.jl`.
