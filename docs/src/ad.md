# Automatic differentiation

Prices are generic in the number type of the inputs: the model, the underlying, the curve and the maturity. Passing an
automatic differentiation number gives the derivatives of the price, of any order, with each back-end loaded with its
package:

| Method | DualNumbers | HyperDualNumbers | TaylorSeries | ForwardDiff | ReverseDiff | Mooncake |
|:--|:--|:--|:--|:--|:--|:--|
| `Vibrato`, Gaussian last step | yes | yes | yes | yes | yes | reverse |
| `Vibrato`, FFT density | yes | yes | yes | yes | yes | forward |
| `VibratoSaltando` | yes | yes | yes | yes | yes | reverse |
| `PureLRM` | yes | yes | yes | yes | yes | not tested |
| American options | yes | yes | yes | yes | yes | forward on Lévy processes |

"Reverse" and "forward" are the modes the tests check for Mooncake; the other entries are tested against references or
against each other.

ForwardDiff covers nested duals, so `ForwardDiff.hessian` gives the second derivatives of several inputs at once:

```@example ad
using VibratoMonteCarlo, FinancialMonteCarlo, ForwardDiff

mc = MonteCarloConfiguration(20_000, 10; sampling = FinancialMonteCarlo.AntitheticMC())
price(x) = pricer(BlackScholesProcess(x[2], Underlying(x[1], 0.0)), ZeroRate(0.02), Vibrato(mc, VibratoMonteCarloGaussLegendre()), BinaryEuropeanOption(1.0, 100.0))
ForwardDiff.hessian(price, [100.0, 0.2])
```

TaylorSeries gives every derivative up to an order in one evaluation, `Taylor1` for one input, `TaylorN` for several:

```@example ad
using TaylorSeries

S0 = Taylor1([100.0, 1.0], 4)   # S0 + t, up to t^4
p = pricer(BlackScholesProcess(0.2, Underlying(S0, 0.0)), ZeroRate(0.02), PureLRM(mc, VibratoMonteCarloGaussLegendre()), EuropeanOption(1.0, 100.0))
[factorial(k) * p[k] for k = 0:4]   # value, delta, gamma, speed, fourth derivative
```

## How it works

The sample points of the inner integration are built from the values of the inputs, so the payoff is never
differentiated, and every derivative enters through the density, in the factor ``m(g) = e^{g - \mathrm{v}(g)}``
evaluated at the log density ``g``: it has value one and, by Faà di Bruno's formula, its derivatives are the likelihood
ratio weights of every order. An extension implements the value ``\mathrm{v}`` for each number type, see
[Why automatic differentiation works](@ref).

Mooncake differentiates the plain `Float64` call instead: its extension makes ``\mathrm{v}`` a constant for
differentiation and gives ``m`` its exponential form. Reverse mode works on the Gaussian methods and on `VibratoSaltando`
(with ``\kappa`` and ``T``); the FFT `Vibrato` works in forward mode, with `Mooncake.Config(enable_nfwd = false)`. The
first call compiles for minutes, so for few inputs ForwardDiff is the better choice.

## GPU

With CUDA.jl loaded, a `MonteCarloConfiguration` with `FinancialMonteCarlo.CudaMode()` simulates the paths and integrates
the last step on the GPU:

```julia
using CUDA
mc = MonteCarloConfiguration(100_000, 10; sampling = FinancialMonteCarlo.StandardMC(), mode = FinancialMonteCarlo.CudaMode())
pricer(model, zero_rate, Vibrato(mc, VibratoMonteCarloGaussLegendre()), option)
```

- On the device: the paths of FinancialMonteCarlo.jl, the variance of Heston, the subordinator ``G_T`` of VG and NIG with
  plain Monte Carlo (FinancialMonteCarlo.jl's Gamma and inverse Gaussian sampler), the paths of Kou and Merton up to their
  last jump, the FFT density of the last step (`Vibrato` and `PureLRM` on Lévy processes, computed by FinancialFFT.jl in
  the mode of the configuration), and the inner integration on the matrix of paths and nodes, with the samples of
  `VibratoMonteCarloStandard` drawn by CURAND.
- On the host, then copied to the device: the antithetic ``G_T`` (pairs of quantiles) and the jump times and sizes of Kou
  and Merton (drawn by FinancialMonteCarlo.jl, as for its own GPU paths).
- On the host: American options and `PureLRM` on Itô processes.
- The numbers must be isbits on the device: hyper duals and ForwardDiff (nested too) are tested. TaylorSeries and
  ReverseDiff cannot run there.
