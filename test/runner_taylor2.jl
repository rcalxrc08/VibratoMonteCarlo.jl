using TaylorSeries, FinancialMonteCarlo, FinancialToolbox, VibratoMonteCarlo, HyperDualNumbers, FinancialFFT, FastGaussQuadrature
#S0 = 100.0;
S0 = taylor_expand(identity, 100.0, order = 5);
K = 105.0;
r = 0.02;
T = 1.0;
d = 0.0;
D = 90.0;

Nsim = 10_000;
Nstep = 300;
sigma = 0.2;
p = 0.3;
lam = 50.0;
lamp = 30.0;
lamm = 20.0;
mu1 = 0.03;
sigma1 = 0.02;
mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
# mc = MonteCarloConfiguration(Nsim, Nstep);

rfCurve = ZeroRate(r);

EUData = EuropeanOption(T, K)
Model = KouProcess(sigma, lam, p, lamp, lamm, Underlying(S0, d));
A = 600.0;
N = 20000;

method = LewisMethod(A, N);

@show result1 = vibrato_saltando(Model, rfCurve, mc, EUData, VibratoMonteCarlo.VibratoMonteCarloGaussHermite(200));
@show result4 = [derivative(pricer(Model, rfCurve, method, EUData), i)[0] for i = 0:5];
@show result5 = pricer(Model, rfCurve, mc, EUData);
