# Benchmarks

A timing script run by hand, not by CI nor the tests.

## Inputs

| Load | Why |
|:--|:--|
| `bench.jl` | `@btime` of `Vibrato`, `VibratoSaltando` and FinancialMonteCarlo's own pricer on Kou and Black–Scholes, hyper dual `S0` |
| `Project.toml` | its environment: the package is added with `Pkg.develop(path = pwd())`, its compat follows the package's |
| `references/workflows.md` | the command, compat bumps in every sub-environment |
| `references/gpu.md` | the GPU timings live in `test/cuda/timings.jl`, not here |

## Process

Run one benchmark per fresh process; say on which machine and Julia version a number was taken when it goes into a
reference or the docs. A FinancialMonteCarlo compat bump goes into `Project.toml` here too. The script still prices an
array of options, which the package no longer supports (#72).

## Outputs

Timings printed to the terminal; nothing is written to the repository.
