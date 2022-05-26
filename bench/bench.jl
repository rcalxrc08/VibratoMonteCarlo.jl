using DualNumbers, HyperDualNumbers, FinancialMonteCarlo, Statistics, BenchmarkTools
include("vibrato_saltando.jl");
include("vibrato.jl");
@show "KouModel"
# S0 = 100.0;
# S0 = dual(100.0, 1.0);
S0 = hyper(100.0, 1.0, 1.0, 0.0);
K = 100.0;
r = 0.02;
# r = hyper(0.02, 1.0, 1.0, 0.0);
T = 1.0;
d = 0.01;
D = 90.0;

Nsim = 10000;
Nsim2 = 100;
Nstep = 30;
# sigma = dual(0.2, 1.0);
sigma = 0.2;
p = 0.3;
lam = 5.0;
# lam = dual(5.0, 1.0);
lamp = 30.0;
# lamp = dual(30.0, 1.0);
lamm = 20.0;
mc = MonteCarloConfiguration(Nsim, Nstep);
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());

rfCurve = ZeroRate(r);

EUData = BinaryEuropeanOption(T, K)
Model = KouProcess(sigma, lam, p, lamp, lamm, Underlying(S0, d));
@btime result = vibrato_saltando(Model, rfCurve, mc, EUData, VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0));
@btime EuPrice = pricer(Model, rfCurve, mc, EUData);

bs = BlackScholesProcess(sigma, Underlying(S0, d));
@btime result = vibrato(bs, rfCurve, mc, EUData, VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0));
@btime EuPrice = pricer(bs, rfCurve, mc, EUData);
