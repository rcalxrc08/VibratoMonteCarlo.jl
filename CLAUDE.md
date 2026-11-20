# CLAUDE.md

VibratoMonteCarlo.jl computes option prices **and sensitivities of any order** with Monte Carlo methods that stay smooth
(C^∞) under automatic differentiation: AD numbers (DualNumbers, HyperDualNumbers, ForwardDiff, ReverseDiff,
TaylorSeries) go in as model or market inputs and the derivatives come out of the price; Mooncake differentiates the
plain `Float64` call. The methods plug into FinancialMonteCarlo's `pricer(model, zero_rate, method, payoff)`, as
FinancialFFT's `LewisMethod` does; paths come from FinancialMonteCarlo, densities from FinancialFFT, closed forms from
FinancialToolbox (`references/dependencies.md`, with the private registry to add before `Pkg.instantiate()`).

This workspace follows the Interpretable Context Methodology (layered context, arXiv:2603.16021): load only what the
task needs.

| Layer | Where | What |
|:--|:--|:--|
| 0 | this file | identity, map, rules that hold everywhere |
| 1 | `CONTEXT.md` | which stage or folder a task belongs to, shared references |
| 2 | `src/0N_*/CONTEXT.md`, `ext/`, `test/`, `docs/` | a contract per stage: inputs to read, process, outputs and tests |
| 3 | `references/*.md` | stable facts and measurements, loaded when a contract lists them |

Map: `src/VibratoMonteCarlo.jl` includes the stages in order, `01_methods` → `02_simulation` → `03_inner_integration` →
`04_european` → `05_american`; `ext/` holds one extension per AD backend; `test/` the suite (`test/cuda/` the opt-in GPU
tests); `docs/` the Documenter site; `bench/` a benchmark script; `examples/` the README examples; `etc/` the logo and
`readme.jl`. `src/deps/` and other untracked files are experiments, not package code (`references/workflows.md`).

Rules for every change (details in `references/conventions.md` and `references/autodiff.md`):

- Sample points and payoffs from values only (`v_value`); every derivative flows through the density weight `v_mod`.
- Generic number types: arrays built from the type of the inputs, no Float64 constants on the pricing path; Float32 in,
  Float32 out (#45, #70).
- Every pricer calls `set_seed!` first and stays inferred; validation by `ArgumentError`/`DomainError`, never `@assert`.
- Reuse the author's packages (FinancialMonteCarlo for payoffs, paths and laws; FinancialFFT for densities and drifts;
  FinancialToolbox's `blsprice`/`blsbin`; PatternVectors for quadrature weights) and established libraries;
  `@muladd` for every `a * b + c`.
- Format with the repo's `.JuliaFormatter.toml`; update the CONTEXT.md and references touched by a change in the same
  commit; push to `origin` (GitHub) only; commit messages end with the session's `Co-Authored-By` line.
