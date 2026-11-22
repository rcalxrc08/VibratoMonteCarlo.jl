# Decisions

_Architecture decisions, oldest first. Append; a reversed decision gets a new entry that says "supersedes D-NNN"._
(Layer 3 reference.) Read before a refactor: a structure that looks accidental may be an entry here. Details stay in
the reference each entry points to.

## D-001 (by 2026-10-28, 86b4fc6) Sample from values, differentiate through the density weight

- **Decision:** sample points and payoffs see only `v_value` (every derivative level stripped); every derivative flows
  through the density weight `v_mod`.
- **Why:** the estimator must stay smooth in the inputs: derivatives of the sampled points would bring back the
  discontinuities of the payoff.
- **Consequences:** a new AD backend is an extension adding `v_value`/`v_mod` methods (`autodiff.md`); rule in
  `CLAUDE.md`.
- **Status:** active

## D-002 (2026-11-01, f53f13b; 2026-11-21, 4b8d89f) Smiles are priced strike by strike, no arrays of options

- **Decision:** pricers dispatch on smiles directly (`SmileStrike{O}(K, isCall)`, `map_strikes`), with no intermediate
  vector of option objects.
- **Why:** smiles 200 times lighter (#42 port), and the inner integration shares the draws of a path across strikes.
- **Consequences:** `european.md`, `public_api.md`.
- **Status:** active

## D-003 (2026-11-19, a25cdcc; #58, 0.6.0) Explicit exceptions, never `@assert`

- **Decision:** `check_argument` (`ArgumentError`) in method constructors, `check_domain` (`DomainError`) for values out
  of their domain.
- **Why:** `@assert` throws a generic error and may be disabled; FinancialCalibrator turns every `DomainError` of the
  loss into `Inf`, so a misconfigured method must fail before any pricing.
- **Consequences:** `conventions.md`.
- **Status:** active

## D-004 Every pricer reseeds first and stays inferred

- **Decision:** each public pricer calls `FinancialMonteCarlo.set_seed!(mcBaseData)` before drawing; pricers are
  `@inferred` in the tests.
- **Why:** reproducible prices, and derivatives computed by different AD passes on the same paths.
- **Consequences:** `conventions.md`, `testing.md`.
- **Status:** active

## D-005 (2026-11-20, 9d24eac; merged in PR #73 on 2026-11-21) Layered context (Interpretable Context Methodology)

- **Decision:** `CLAUDE.md` (L0), `CONTEXT.md` routing (L1), stage contracts (L2), `references/` (L3), `state/` (L4,
  added 2026-11-21); `test/context.jl` checks the module's include list against the contracts.
- **Why:** same layout as FinancialCalibrator; one contract per stage loads only what a task needs.
- **Consequences:** a change updates the contracts, references and `state/STATUS.md` it touches in the same commit.
- **Status:** active

## D-006 (2026-11-21, 698ddbd; #74) Weights that do not depend on the path are computed once per price

- **Decision:** the Analytic rule precomputes nodes and `trapezoid weights × φ(z) × step/2` once per price.
- **Why:** 4b8d89f cut allocations 20000 → 24 but made the Analytic binary 2x slower (an `exp` per node and path);
  allocation counts are not a proxy for speed.
- **Consequences:** hot-loop changes are timed before and after in fresh processes (`bench/CONTEXT.md`).
- **Status:** active

## D-007 (2026-11-22, #68) Rates and dividends through their integrals

- **Decision:** every drift over an interval and every discount goes through `interval_integral`, `add_drift`,
  `discount_factor` and `mean_rate` (`src/01_methods/utils.jl`), exact for numbers and FinancialMonteCarlo curves; with
  numbers they are the constant-rate formulas (dispatch), so constant-rate prices do not change. Code that needs one
  rate for every step (Heston's Euler scheme, the frozen laws of American options) refuses curves with `constant_rate`.
- **Why:** FinancialMonteCarlo 0.9 rates and dividends can be curves; `r * T` and `exp(-r * T)` on a curve failed with a
  `MethodError` in every European method.
- **Consequences:** rule in `CLAUDE.md`; `references/conventions.md`; `test/rate_curves.jl`.
- **Status:** active
