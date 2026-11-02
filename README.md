# VibratoMonteCarlo.jl <img src="etc/logo.png" width="40">  

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rcalxrc08.github.io/VibratoMonteCarlo.jl/)
[![CI](https://github.com/rcalxrc08/VibratoMonteCarlo.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/rcalxrc08/VibratoMonteCarlo.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/rcalxrc08/VibratoMonteCarlo.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rcalxrc08/VibratoMonteCarlo.jl)
##### This is a Julia package containing some useful Financial function for Pricing and Risk Management for Equity products.

This package is built on top of [FinancialMonteCarlo.jl](https://github.com/rcalxrc08/FinancialMonteCarlo.jl) and [FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl).
Standard montecarlo methods lacks of differentiability, which makes automatic differentiation useless.
The main aim of this package is to provide a feasible way to compute sensitivities of any order for various types of payoffs using montecarlo methods.
It currently contains the following capabilities:

- Support for the following Single Name Models:
    - Black Scholes
    - Kou
    - Merton
    - Normal Inverse Gaussian
    - Variance Gamma
- Support for the following payoffs:
    - European Options
    - Binary European Options
	- American Options (WIP)
- Support for the following metrics:
    - likelihood ratio method
    - vibrato montecarlo method (for any type of process, not just Ito).
	- vibrato saltando montecarlo method (for Finite Activity Levy).
- Support for the following integration methods:
    - pure numerical
    - montecarlo (standard, antithetic, sobol where feasible).
    - gaussian hermite (where feasible)
	
- Currently supports:
	- [DualNumbers.jl](https://github.com/JuliaDiff/DualNumbers.jl)
	- [HyperDualNumbers.jl](https://github.com/JuliaDiff/HyperDualNumbers.jl)
	- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)
	- [ReverseDiff.jl](https://github.com/JuliaDiff/ReverseDiff.jl)
	- [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl)

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
[examples/hello_world.jl](examples/hello_world.jl) prices a european call option with underlying varying according to the Black Scholes Model, given the volatility, and its delta. With FinancialMonteCarlo and DualNumbers installed, run it with:
```julia
using VibratoMonteCarlo
include(joinpath(pkgdir(VibratoMonteCarlo), "examples", "hello_world.jl"))
```