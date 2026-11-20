# American options

_Stage 05 reference._ (Layer 3 reference: stable facts, loaded only when a stage CONTEXT.md lists it.)

`src/05_american/vibrato_am.jl`, `pricer(process, rfCurve, Vibrato(mc, vb_mc), am_opt::AmericanPayoff)`: exercise at the
`Nstep` grid dates, policy from Longstaff–Schwartz (regression on `1, S, S²` of the in-the-money paths, `S` scaled by
`S0`, `v_value`s only). Two methods.

## Black–Scholes: vibrato over every exercise date

- Pass 1 (independent paths, values only, `exercise_regions`) turns the regressions into exercise regions in log price
  (`positive_intervals`), fixed in price space.
- Pass 2 simulates `2 Nstep` half steps (same MC configuration, shared rng) so that every exercise date `t_k` is the
  midpoint of `[s_k, s_{k+1}]`, `s_k = (k - 1/2) dt`, where `log S_{t_k}` is a Brownian bridge (mean of the ends, std
  `σ √dt / 2`). The configuration is rebuilt with the field constructor
  `MonteCarloConfiguration(Nsim, 2Nstep, sampling, mode, rng, reseed)`, which keeps the random state; if
  FinancialMonteCarlo ever hides it, use `Accessors.@set`.
- Backwards over the dates for all the paths at once (`bridge_step`: `bridge_vibrato` path by path on the host, a
  broadcast per piece on the device), `bridge_vibrato` integrates "exercise if in the region, else continue with the path
  value" against the bridge density (`lrm_segment`: Gauss–Legendre on each region piece, points from values,
  derivatives in `v_mod`); the payoff at `T` gets the usual vibrato over the last half step.
- Since the optimal boundary does not depend on `S0`, a policy fixed in price space gives the `S0` derivatives of every
  order (the regression error enters them only at second order).
- On Julia 1.10 a generator whose elements build a closure is not inferred, so the bridge sums its segments with
  `sum(f, segments)` over callable structs (`PayoffAt`, `UnitPayoff`), not
  `sum(lrm_segment(…, x -> payout_f(x, am_opt)) for …)`.

## Lévy and Heston: frozen decisions

The decision to exercise at `t_k` is taken at `t_{k-1}` and the cash flow is the vibrato expectation
`E[payoff(S_{t_k}) | S_{t_{k-1}}]` over one step (`frozen_laws` returns the paths and the laws of one step: `FrozenFFT`,
the FFT `FFTStepLaw` and its `value_law` for the values only regressions, which needs `VibratoMonteCarloAnalytic`, for
Lévy; `FrozenHeston`, Gaussian with the variance of the path from `simulate_with_vol`, for Heston; `step_payoffs`
smooths a set of paths at once). One AD cash flow per path at its exercise date. Other `ItoProcess`es throw.

## Limits

The gamma of the frozen method misses the exercise boundary term; the BS bridge gamma is unbiased but its variance grows
with `Nstep`: numbers in `references/method_limits.md`. Never add a smoothing of the decisions to fix it. American
pricers return Float64 for Float32 models. GPU split (regressions on the host, the rest on the device):
`references/gpu.md`.
