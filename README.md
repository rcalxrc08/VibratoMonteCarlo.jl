# VibratoMonteCarlo.jl <img src="etc/logo.png" width="40">  

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rcalxrc08.github.io/VibratoMonteCarlo.jl/)
[![CI](https://github.com/rcalxrc08/VibratoMonteCarlo.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/rcalxrc08/VibratoMonteCarlo.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/rcalxrc08/VibratoMonteCarlo.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rcalxrc08/VibratoMonteCarlo.jl)
##### This is a Julia package containing some useful Financial function for Pricing and Risk Management for Equity products.

This package is built on top of [FinancialMonteCarlo.jl](https://github.com/rcalxrc08/FinancialMonteCarlo.jl) and [FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl).
Standard montecarlo methods lacks of differentiability, which makes automatic differentiation useless.
The main aim of this package is to provide a feasible way to compute sensitivities of any order for various types of payoffs using montecarlo methods.
The [documentation](https://rcalxrc08.github.io/VibratoMonteCarlo.jl/) has the theory of the methods, examples and the API.
It currently contains the following capabilities:

- Support for the following Single Name Models:
    - Black Scholes
    - Heston
    - Kou
    - Merton
    - Normal Inverse Gaussian
    - Variance Gamma
- Support for the following payoffs:
    - European Options
    - Binary European Options
    - Smiles of European and Binary European Options (one maturity, many strikes)
    - American Options (exercise at the time steps, Longstaff-Schwartz policy)
- Support for the following pricing methods, passed to FinancialMonteCarlo's `pricer(model, zero_rate, method, payoff)` (as FinancialFFT's `LewisMethod`):
    - `PureLRM(mc, vb_mc)`: likelihood ratio method on the density at maturity (Gaussian or from the FFT).
    - `Vibrato(mc, vb_mc)`: vibrato montecarlo method, for Black Scholes, Heston and any Levy process (FFT density of the last step).
    - `VibratoSaltando(mc, vb_mc)`: vibrato conditioned on the last jump time (Kou, Merton) or on the subordinator (Variance Gamma, Normal Inverse Gaussian, with derivatives over `κ` and `T` too); exact for any number of steps.
    - `Vibrato` on American options: vibrato over every exercise date in Black Scholes, frozen exercise decisions for Levy and Heston.
- Support for the following integration methods of the last step:
    - `VibratoMonteCarloGaussLegendre`: Gauss-Legendre split at the strikes, the most accurate.
    - `VibratoMonteCarloAnalytic`: quadrature on a uniform grid.
    - `VibratoMonteCarloStandard`: montecarlo (standard, antithetic, sobol).
- Support for the following automatic differentiation backends:
	- [DualNumbers.jl](https://github.com/JuliaDiff/DualNumbers.jl)
	- [HyperDualNumbers.jl](https://github.com/JuliaDiff/HyperDualNumbers.jl)
	- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) (nested too, e.g. `ForwardDiff.hessian`)
	- [ReverseDiff.jl](https://github.com/JuliaDiff/ReverseDiff.jl)
	- [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl)
	- [Mooncake.jl](https://github.com/chalk-lab/Mooncake.jl) (on the plain `Float64` call)
- GPU support with `CudaMode` of FinancialMonteCarlo (`using CUDA`), for hyper duals and ForwardDiff.

To be noticed that the methods are $C^{\infty}$, not as standard montecarlo methods.
## How to Install
To install the package simply type on the Julia REPL the following:
```julia
] add VibratoMonteCarlo
```
## How to Test
After the installation, to test the package type on the Julia REPL the following:
```julia
] test VibratoMonteCarlo
```
## Hello World: Pricing European Call Option in Black Scholes Model
The following example ([examples/hello_world.jl](examples/hello_world.jl)) shows how to price a european call option with underlying varying according to the Black Scholes Model, given the volatility, and its delta:
<!-- BEGIN examples/hello_world.jl -->
```julia
# Hello World: pricing a European call option in the Black Scholes model, with its delta.
#Import the Package
using VibratoMonteCarlo, FinancialMonteCarlo, DualNumbers;

#Define Spot Datas
S0 = dual(100.0, 1.0);
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;
#Define VibratoMonteCarlo Parameters
Nsim = 10000;
Nstep = 30;
#Define Model Parameters
σ = 0.2;
#Build the Structs
A = 600.0;
N = 18;
mc = MonteCarloConfiguration(Nsim, Nstep);
zeroRate = ZeroRate(r);
underlying = Underlying(S0, d); #Underlying relative data

#Define The Option
EU_payoff = EuropeanOption(T, K)
#Define the Model
Model = BlackScholesProcess(σ, underlying);

#Price and sensitivity
@show EuPrice = pricer(Model, zeroRate, PureLRM(mc, VibratoMonteCarlo.VibratoMonteCarloAnalytic(20000, -9.0, 19.0)), EU_payoff);
```
<!-- END examples/hello_world.jl -->

With FinancialMonteCarlo and DualNumbers installed, run it with:
```julia
using VibratoMonteCarlo
include(joinpath(pkgdir(VibratoMonteCarlo), "examples", "hello_world.jl"))
```