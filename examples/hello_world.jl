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
@show EuPrice = pure_lrm(Model, zeroRate, mc, EU_payoff, VibratoMonteCarlo.VibratoMonteCarloAnalytic(20000, -9.0, 19.0));
