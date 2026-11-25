# CLAUDE.md

VibratoMonteCarlo.jl computes option prices **and sensitivities of any order** with Monte Carlo methods that stay smooth
(C^∞) under automatic differentiation: AD numbers (DualNumbers, HyperDualNumbers, ForwardDiff, ReverseDiff,
TaylorSeries) go in as model or market inputs and the derivatives come out of the price; Mooncake differentiates the
plain `Float64` call. The methods plug into FinancialMonteCarlo's `pricer(model, underlying, zero_rate, method, payoff)`, as
FinancialFFT's `LewisMethod` does; paths come from FinancialMonteCarlo, densities from FinancialFFT, closed forms from
FinancialToolbox (`references/dependencies.md`, with the private registry to add before `Pkg.instantiate()`). Goals, success criteria
and non-goals: `references/mission.md`; why things are the way they are: `references/decisions.md`.

This workspace follows the Interpretable Context Methodology (layered context, arXiv:2603.16021): load only what the
task needs.

| Layer | Where | What |
|:--|:--|:--|
| 0 | this file | identity, map, rules that hold everywhere |
| 1 | `CONTEXT.md` | which stage or folder a task belongs to, shared references |
| 2 | `src/0N_*/CONTEXT.md`, `ext/`, `test/`, `docs/`, `etc/`, `bench/` | a contract per stage: inputs to read, process, outputs, runtime checks |
| 3 | `references/*.md` | stable facts, decisions and measurements, loaded when a contract lists them |
| 4 | `state/STATUS.md` | what is in flight: current focus, open issues by priority, known problems |

Startup: read `state/STATUS.md`, then route the task through `CONTEXT.md`. Runtime truth comes from a Kaimon session
for this project (`kaimon.toml` pins the Julia); each contract's Verify table lists the checks to run after a change.

Map: `src/VibratoMonteCarlo.jl` includes the stages in order, `01_methods` → `02_simulation` → `03_inner_integration` →
`04_european` → `05_american`; `ext/` holds one extension per AD backend; `test/` the suite (`test/cuda/` the opt-in GPU
tests); `docs/` the Documenter site; `bench/` a benchmark script run by hand; `examples/` the README examples, copied into
`README.md` by `etc/readme.jl`; `etc/` also the logo. `src/deps/` and other untracked files are experiments, not package code (`references/workflows.md`).

Rules for every change (details in `references/conventions.md` and `references/autodiff.md`):

- Sample points and payoffs from values only (`v_value`); every derivative flows through the density weight `v_mod`.
- Generic number types: arrays built from the type of the inputs, no Float64 constants on the pricing path; Float32 in,
  Float32 out (#45, #70).
- Rates and dividends are numbers or curves: drifts and discounts through `interval_integral`/`add_drift`/
  `discount_factor` (`src/01_methods/utils.jl`, #68), never `r * T` on a rate.
- Every pricer calls `set_seed!` first and stays inferred; validation by `ArgumentError`/`DomainError`, never `@assert`.
- Reuse the author's packages (FinancialMonteCarlo for payoffs, paths and laws; FinancialFFT for densities and drifts;
  FinancialToolbox's `blsprice`/`blsbin`; PatternVectors for quadrature weights) and established libraries;
  `@muladd` for every `a * b + c`.
- Format with the repo's `.JuliaFormatter.toml`; update the CONTEXT.md, references and `state/STATUS.md` touched by a
  change in the same commit (a new rule also gets an entry in `references/decisions.md`); push to `origin` (GitHub) only; commit messages end with the session's `Co-Authored-By` line.
