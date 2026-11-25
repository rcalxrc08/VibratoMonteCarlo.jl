# Mission

_Why the package exists and how to tell whether it succeeds._ (Layer 3 reference: stable, changes only with an entry in
`decisions.md`.)

## Purpose

Monte Carlo prices are not differentiable for discontinuous payoffs (binaries, barriers) and give noisy higher-order
derivatives. VibratoMonteCarlo computes option prices **and sensitivities of any order** with Monte Carlo methods that
stay smooth (C^∞) under automatic differentiation: the vibrato step integrates the last step analytically or by
quadrature, so AD numbers that go in as model or market inputs come out as exact derivatives of the estimator.

## Goals

- Sensitivities of any order (DualNumbers, HyperDualNumbers, ForwardDiff, ReverseDiff, TaylorSeries) of the price of
  European and American options, for Black–Scholes, Heston and the Lévy models of FinancialMonteCarlo.
- Mooncake differentiates the plain `Float64` call.
- The methods plug into FinancialMonteCarlo's `pricer(model, underlying, zero_rate, method, payoff)`, so FinancialCalibrator can
  calibrate with them.
- CPU and CUDA from the same code; Float32 in, Float32 out.

## Constraints

- Paths, payoffs and laws from FinancialMonteCarlo; densities and drifts from FinancialFFT; closed forms from
  FinancialToolbox; quadrature weights from PatternVectors (`dependencies.md`).
- Every pricer reseeds first (`set_seed!`) and is inferred (D-004).
- Julia ≥ 1.10 (`Project.toml`); CI on Julia 1.

## Success criteria

- `test/` reproduces closed-form prices and Greeks (Black–Scholes) and finite-difference references within the
  tolerances of `references/testing.md`; the CUDA tests (`test/cuda/`) agree with the host.
- `references/method_limits.md` states, per model and method, which derivatives are supported and why the others are
  not.
- Hot-path changes carry before/after timings (`bench/`).

## Non-goals

- Simulation of models and payoffs: FinancialMonteCarlo. Characteristic functions and Fourier pricing: FinancialFFT.
- Calibration: FinancialCalibrator.
- Unbiased derivatives of discontinuous payoffs by plain pathwise AD: that is exactly what the vibrato step replaces.
