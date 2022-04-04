using DualNumbers, HyperDualNumbers, FinancialMonteCarlo, FinancialFFT, Statistics, Distributions, FinancialToolbox
include("vibrato_saltando.jl");
include("vibrato.jl");
@show "KouModel"
# S0 = 100.0;
#S0 = dual(100.0, 1.0);
S0 = hyper(100.0, 1.0, 1.0, 0.0);
K = 80.0;
r = 0.01;
# r = dual(0.02, 1.0);
# r = hyper(0.02, 1.0, 1.0, 0.0);
T = 1.0;
#T = dual(1.0, 2.0);
d = 0.01;
D = 90.0;

Nsim = 100_000;
Nsim2 = 200;
Nstep = 300;
# sigma = hyper(0.2, 1.0, 1.0, 0.0);
# sigma = dual(0.2, 1.0);
sigma = 0.2;
p = 0.3;
lam = 5.0;
# lam = dual(5.0, 1.0);
lamp = 30.0;
# lamp = dual(30.0, 1.0);
# lamp = hyper(30.0, 1.0, 1.0, 0.0);
lamm = 20.0;
mu1 = 0.03;
sigma1 = 0.02;
mc = MonteCarloConfiguration(Nsim, Nstep);
#mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
# mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SobolMode());
method_lewis = LewisMethod(400.0, 200000);
rfCurve = ZeroRate(r);

EUData = EuropeanOption(T, K)
Model = KouProcess(sigma, lam, p, lamp, lamm, Underlying(S0, d));
# Model = MertonProcess(sigma, lam, mu1, sigma1, Underlying(S0, d));
@show result = vibrato_saltando(Model, rfCurve, mc, EUData, VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0));
@show EuPrice = pricer(Model, rfCurve, mc, EUData);
@show EuPrice2 = pricer(Model, rfCurve, method_lewis, EUData);

bs = BlackScholesProcess(sigma, Underlying(S0, d));
@show result = vibrato(bs, rfCurve, mc, EUData, VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0));
@show EuPrice = pricer(bs, rfCurve, mc, EUData);
@show EuPrice2 = pricer(bs, rfCurve, method_lewis, EUData);
@show blsprice(S0, K, r, T, sigma, d);

sigma_zero = 0.2;
kappa = 0.2;
theta = 0.2;
lambda = 0.0;
rho_ = 0.2;
hest = HestonProcess(sigma, sigma_zero, lambda, kappa, rho_, theta, Underlying(S0, d));
@show result = vibrato(hest, rfCurve, mc, EUData, VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0));
@show EuPrice = pricer(hest, rfCurve, mc, EUData);


vg = VarianceGammaProcess(sigma, 0.1, 0.01, Underlying(S0, d));
@show result = vibrato_saltando(vg, rfCurve, mc, EUData, VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0));
@show EuPrice = pricer(vg, rfCurve, mc, EUData);
@show EuPrice2 = pricer(vg, rfCurve, method_lewis, EUData);

nothing