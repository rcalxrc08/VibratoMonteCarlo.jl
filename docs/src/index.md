```@raw html
<img class="pkg-hero-logo" src="assets/logo.png" alt="VibratoMonteCarlo.jl logo" width="96" height="96">
```

# VibratoMonteCarlo.jl

```@raw html
<p class="pkg-roles"><b>Monte Carlo Greeks</b> · <b>Vibrato</b> · <b>Automatic differentiation</b> · <b>GPU</b></p>
<div class="pkg-intro">
<p>Plain Monte Carlo estimators are not differentiable: the pathwise derivative of a binary option is zero, and
automatic differentiation gives nothing useful on discontinuous payoffs. VibratoMonteCarlo.jl builds estimators that are
smooth in every parameter, likelihood ratio, Vibrato and Vibrato Saltando, so that sensitivities of any order come out of
automatic differentiation of the price, for Black and Scholes, Heston and Lévy models, European, binary and American
options.</p>
</div>
<nav class="pkg-links">
<a class="pkg-btn" href="https://github.com/rcalxrc08/VibratoMonteCarlo.jl"><svg viewBox="0 0 16 16" aria-hidden="true"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"/></svg> GitHub</a>
<a class="pkg-btn" href="https://github.com/rcalxrc08/VibratoMonteCarlo.jl/stargazers">★ Star</a>
<a class="pkg-btn" href="#Quick-start">Get started ↓</a>
<a class="pkg-btn" href="https://rcalxrc08.github.io/">All projects</a>
</nav>
```

![Delta and gamma of a binary call with Vibrato against Black and Scholes](assets/figures/binary_greeks-light.svg)![Delta and gamma of a binary call with Vibrato against Black and Scholes](assets/figures/binary_greeks-dark.svg)

## Relation to FinancialMonteCarlo.jl

The package adds methods to the `pricer` of [FinancialMonteCarlo.jl](https://github.com/rcalxrc08/FinancialMonteCarlo.jl),
as [FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl) does, and uses its types for everything else: the
model (`BlackScholesProcess`, `HestonProcess`, …), the underlying (`Underlying`), the curve (`ZeroRate`), the payoff
(`EuropeanOption`, `BinaryEuropeanOption`, `AmericanOption`, and the smiles) and the `MonteCarloConfiguration`. The
method wraps the configuration, so moving from plain Monte Carlo to Vibrato changes one argument:

```julia
pricer(model, zero_rate, mc, option)                                            # Monte Carlo, FinancialMonteCarlo.jl
pricer(model, zero_rate, Vibrato(mc, VibratoMonteCarloGaussLegendre()), option)  # Vibrato, VibratoMonteCarlo.jl
pricer(model, zero_rate, LewisMethod(400.0, 200_000), option)                   # Fourier, FinancialFFT.jl
```

The implementation is based on the Master Thesis of Nicola Scaramuzzino [scaramuzzino2017](@cite) and the subsequent
development.

## Features

- Pricing methods: [`Vibrato`](@ref), [`VibratoSaltando`](@ref) and [`PureLRM`](@ref), with American options through `Vibrato`, see [Methods](methods.md). The theory is in [Pathwise and likelihood ratio methods](lrm.md), [Conditioning: Vibrato and Vibrato Saltando](vibrato.md), [Numerical refinements](refinements.md) and [American options](american.md).
- Inner integrations: [`VibratoMonteCarloGaussLegendre`](@ref) (split at the strikes, exact to about 1e-14), [`VibratoMonteCarloAnalytic`](@ref) and [`VibratoMonteCarloStandard`](@ref) (Monte Carlo: standard, antithetic, Sobol).
- Models: Black and Scholes, Heston, Kou, Merton, Variance Gamma, Normal Inverse Gaussian. The densities of the Lévy processes come from FinancialFFT.jl.
- Payoffs: European and binary options, calls and puts, their smiles (several strikes, one maturity), and American options.
- Sensitivities of any order with DualNumbers.jl, HyperDualNumbers.jl, ForwardDiff.jl (nested too), ReverseDiff.jl, TaylorSeries.jl and Mooncake.jl, see [Automatic differentiation](@ref).
- GPU execution through CUDA.jl, see [GPU](@ref).

## Installation

Some dependencies are registered in a registry of their own:

```julia
import Pkg
Pkg.Registry.add(Pkg.RegistrySpec(url = "https://github.com/rcalxrc08/LocalRegistry.git"))
Pkg.Registry.add("General")
Pkg.add("VibratoMonteCarlo")
```

## Quick start

Price, delta and gamma of a European call in Black and Scholes, from one evaluation with hyper dual numbers:

```@example quickstart
using VibratoMonteCarlo, FinancialMonteCarlo, FinancialToolbox, HyperDualNumbers

S0 = hyper(100.0, 1.0, 1.0, 0.0)   # the input to differentiate
K, r, d, T, σ = 100.0, 0.02, 0.01, 1.0, 0.2
model = BlackScholesProcess(σ, Underlying(S0, d))
zero_rate = ZeroRate(r)
mc = MonteCarloConfiguration(100_000, 10; sampling = FinancialMonteCarlo.AntitheticMC())

price = pricer(model, zero_rate, Vibrato(mc, VibratoMonteCarloGaussLegendre()), EuropeanOption(T, K))
(value = price.value, delta = price.epsilon1, gamma = price.epsilon12)
```

against the formula of Black and Scholes:

```@example quickstart
exact = blsprice(S0, K, r, T, σ, d)
(value = exact.value, delta = exact.epsilon1, gamma = exact.epsilon12)
```

The same call works for a binary option, whose pathwise gamma is zero, and for the other models:

```@example quickstart
kou = KouProcess(0.2, 1.0, 0.3, 20.0, 10.0, Underlying(S0, d))
price = pricer(kou, zero_rate, VibratoSaltando(mc, VibratoMonteCarloGaussLegendre()), BinaryEuropeanOption(T, K))
(value = price.value, delta = price.epsilon1, gamma = price.epsilon12)
```

## Testing

```julia
] test VibratoMonteCarlo
```

The GPU tests run in their own environment, on a machine with an NVIDIA GPU:

```bash
julia --project=test/cuda -e 'import Pkg; Pkg.instantiate(); include("test/cuda/runtests.jl")'
```

## Contents

```@contents
Pages = ["lrm.md", "vibrato.md", "refinements.md", "american.md", "methods.md", "ad.md", "api.md"]
Depth = 2
```
