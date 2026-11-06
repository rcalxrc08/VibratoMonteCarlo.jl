# Methods

Every estimator is a method struct passed to `pricer`, in the place of the `MonteCarloConfiguration` of
FinancialMonteCarlo.jl. It holds the configuration of the outer simulation and, where there is one, the inner
integration of the last step:

| Method | Constructor | Conditions on | Models |
|:--|:--|:--|:--|
| [`Vibrato`](@ref) | `Vibrato(mc, vb_mc)` | the path up to the last time step | Black and Scholes, Heston, Lévy processes (FFT density) |
| [`VibratoSaltando`](@ref) | `VibratoSaltando(mc, vb_mc)` | the jumps (Kou, Merton) or the subordinator ``G_T`` (VG, NIG) | Kou, Merton, Variance Gamma, Normal Inverse Gaussian |
| [`PureLRM`](@ref) | `PureLRM(mc, vb_mc)` | nothing: the density at maturity | Black and Scholes, Lévy processes (FFT density) |
| [`VibratoMonteCarlo.Conditional`](@ref) | `VibratoMonteCarlo.Conditional(mc)` | as `Vibrato`, closed form inner price | Black and Scholes |
| [`VibratoMonteCarlo.ConditionalSaltando`](@ref) | `VibratoMonteCarlo.ConditionalSaltando(mc)` | as `VibratoSaltando`, closed form inner price | Kou, Merton, Variance Gamma, Normal Inverse Gaussian |

`Vibrato` prices American options too, when the payoff is an `AmericanOption` or a `BinaryAmericanOption`, see
[American options](american.md). The two conditional methods are not exported: they are the exact inner integration limit
of `Vibrato` and `VibratoSaltando` on the same paths, mostly test oracles, and only price European and binary options.

All of them agree with the closed formula of Black and Scholes, here for a binary call and its delta:

```@example methods
using VibratoMonteCarlo, FinancialMonteCarlo, FinancialToolbox, HyperDualNumbers

S0 = hyper(100.0, 1.0, 1.0, 0.0)
K, r, d, T, σ = 105.0, 0.02, 0.01, 1.0, 0.2
model = BlackScholesProcess(σ, Underlying(S0, d))
zero_rate = ZeroRate(r)
option = BinaryEuropeanOption(T, K)
mc = MonteCarloConfiguration(50_000, 10, FinancialMonteCarlo.AntitheticMC())
exact = blsbin(S0, K, r, T, σ, d)

methods = [Vibrato(mc, VibratoMonteCarloGaussLegendre()), PureLRM(mc, VibratoMonteCarloGaussLegendre()), VibratoMonteCarlo.Conditional(mc)]
errors(price) = (value = price.value - exact.value, delta = price.epsilon1 - exact.epsilon1, gamma = price.epsilon12 - exact.epsilon12)
[nameof(typeof(method)) => errors(pricer(model, zero_rate, method, option)) for method in methods]
```

`PureLRM` on Black and Scholes is one integral against the density at maturity: with Gauss-Legendre it is exact up to
rounding, and the Monte Carlo configuration only provides the seed.

## Inner integrations

The inner expectation of every method but the conditional ones is an integral of the payoff against a density, on the
standardized variable of the last step. The integration is the second argument of the method:

| Inner integration | Constructor | Nodes | Laws |
|:--|:--|:--|:--|
| [`VibratoMonteCarloGaussLegendre`](@ref) | `VibratoMonteCarloGaussLegendre(Npoints = 48)` | Gauss-Legendre on ``[-8.5, 8.5 + σ]`` standard deviations, split at the strikes, `Npoints` per piece | Gaussian |
| [`VibratoMonteCarloAnalytic`](@ref) | `VibratoMonteCarloAnalytic(Npoints = 200, x_min = -7.0, x_max = 7.0)` | uniform grid on ``[x_{min}, x_{max}]`` standard deviations | Gaussian and FFT density |
| [`VibratoMonteCarloStandard`](@ref) | `VibratoMonteCarloStandard(Nsim_lrm = 200)` | Monte Carlo samples, new for every path: plain, antithetic or Sobol, as the `MonteCarloConfiguration` | Gaussian |

The payoff has a kink (a call) or a jump (a binary) at the strike. Split at the strike, every piece of the Gauss-Legendre
rule has a smooth integrand and converges exponentially; across the jump, Gauss-Hermite, a uniform grid and Monte Carlo
converge slowly:

![Relative error of the gamma of a binary call against the number of nodes of one inner integral](assets/figures/quadrature-light.svg)![Relative error of the gamma of a binary call against the number of nodes of one inner integral](assets/figures/quadrature-dark.svg)

Gauss-Legendre is therefore the default choice on Gaussian laws. The FFT density of the Lévy processes, through
`Vibrato`, `PureLRM` and American options, needs `VibratoMonteCarloAnalytic`: it comes from `FFTDensityInverter` of
[FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl), evaluated by `FinancialFFT.density_values` at all the
nodes of a path in one call.

## Choosing a method

- On Black and Scholes and Heston, `Vibrato` with `VibratoMonteCarloGaussLegendre`. On Black and Scholes `PureLRM` is
  exact and needs no simulation.
- On Kou and Merton, `VibratoSaltando`: the last step is exact and at least one time step long. It does not differentiate
  over the jump intensity ``\lambda`` or the maturity ``T``: it warns.
- On Variance Gamma and Normal Inverse Gaussian, `VibratoSaltando`: conditioned on the subordinator, it is exact for any
  `Nstep`, and it differentiates over ``\kappa`` and ``T`` too. The FFT `Vibrato` needs a smooth density of the last
  step, which the Variance Gamma density is only for ``dt / \kappa \gtrsim 2.5``:

![Relative error of the gamma of a Variance Gamma call against the length of the last step](assets/figures/variance_gamma-light.svg)![Relative error of the gamma of a Variance Gamma call against the length of the last step](assets/figures/variance_gamma-dark.svg)

## Smiles

A smile of strikes with one maturity is priced on the same paths, with one price per strike:

```@example methods
smile = FinancialMonteCarlo.EuropeanOptionSmile(T, [90.0, 100.0, 110.0], [false, true, true])
prices = pricer(model, zero_rate, Vibrato(mc, VibratoMonteCarloGaussLegendre()), smile)
[(value = p.value, delta = p.epsilon1, gamma = p.epsilon12) for p in prices]
```

## American options

`Vibrato` prices an option exercisable at the `Nstep` dates of the time grid, with the exercise policy of
Longstaff-Schwartz. On Black and Scholes the policy is fixed in price space and every exercise date is integrated with
the likelihood ratio method over a Brownian bridge, so that the gamma is unbiased for the policy:

```@example methods
put = AmericanOption(T, K, false)
price = pricer(model, zero_rate, Vibrato(MonteCarloConfiguration(20_000, 10), VibratoMonteCarloGaussLegendre()), put)
(value = price.value, delta = price.epsilon1, gamma = price.epsilon12)
```

On Lévy processes and Heston the exercise decisions are frozen when differentiating: value and delta are accurate, the
gamma misses the term of the exercise boundary and is biased low, see [American options](american.md).

## Limits

These are limits of the methods, measured against the Lewis formula, not defects:

- With `VibratoMonteCarloStandard` the gammas are noisier than with quadrature: up to about 4% at ``10^5`` paths and
  `Nstep = 300`, where values and deltas are within 2e-3.
- The bounds of `VibratoMonteCarloAnalytic` must reach far enough in the upper tail, where a call payoff is large: at ±5
  standard deviations the volga of a Variance Gamma call was 3-45% off, at ±7 (the default) it is within noise.
- `VibratoSaltando` on Variance Gamma, with ``a = T / \kappa``: the ``n``-th ``S_0`` derivative of a call at
  ``K = S_0 e^{(r - d - \psi) T}`` exists only for ``a > (n - 1) / 2``, and the estimator has finite variance only for
  ``a > n - 1``. This is the smoothness of the true price, not of the estimator.
- Derivatives over ``\kappa`` through the law of ``G_T`` are the noisiest: at ``5 \cdot 10^5`` paths the standard
  deviation of ``\partial_\kappa`` is 3% on Variance Gamma and 1.2% on NIG.
- `Vibrato` on American options on Black and Scholes has an unbiased gamma, but its variance grows with `Nstep`.
