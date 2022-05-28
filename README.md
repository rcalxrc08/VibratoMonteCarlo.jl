# VibratoMonteCarlo.jl <img src="etc/logo.png" width="40">  

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rcalxrc08.gitlab.io/VibratoMonteCarlo.jl/)
[![pipeline status](https://gitlab.com/rcalxrc08/VibratoMonteCarlo.jl/badges/master/pipeline.svg)](https://gitlab.com/rcalxrc08/VibratoMonteCarlo.jl/commits/master) 
[![codecov](https://codecov.io/gl/rcalxrc08/VibratoMonteCarlo.jl/branch/\x6d6173746572/graph/badge.svg?token=7SOI4KWB60)](https://codecov.io/gl/rcalxrc08/VibratoMonteCarlo.jl)
##### This is a Julia package containing some useful Financial function for Pricing and Risk Management for Equity products.

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
- Support for the following metrics:
    - likelihood ratio method
    - vibrato montecarlo method (for any type of process, not just Ito).
	- vibrato saltando montecarlo method (for Finite Activity Levy).
- Support for the following integration methods:
    - pure numerical
    - montecarlo (standard, antithetic, sobol where feasible).
    - gaussian hermite (where feasible)
	
Currently supports:
	- [DualNumbers.jl](https://github.com/JuliaDiff/DualNumbers.jl), 
	- [HyperDualNumbers.jl](https://github.com/JuliaDiff/HyperDualNumbers.jl)
	- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)
	- [ReverseDiff.jl](https://github.com/JuliaDiff/ReverseDiff.jl)
	- [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl)

To be noticed that the methods are C^inf, not as standard montecarlo methods.
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
The following example shows how to price a european call option with underlying varying according to the Black Scholes Model, given the volatility:
```julia
#Import the Package
using VibratoMonteCarlo,FinancialMonteCarlo, DualNumbers;

#Define Spot Datas
S0=dual(100.0,1.0);
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;
#Define VibratoMonteCarlo Parameters
Nsim=10000;
Nstep=30;
#Define Model Parameters
σ=0.2;
#Build the Structs
A = 600.0;
N = 18;
mc = MonteCarloConfiguration(Nsim, Nstep);
zeroRate=ZeroRate(r);
underlying=Underlying(S0,d); #Underlying relative data

#Define The Option
EU_payoff=EuropeanOption(T,K)
#Define the Model
Model=BlackScholesProcess(σ,underlying);

#Price
@show EuPrice=pure_lrm(Model,zeroRate,mc,EU_payoff, VibratoMonteCarlo.VibratoMonteCarloAnalytic(20000,-9.0,19.0));
```