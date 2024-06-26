using TaylorSeries, FinancialMonteCarlo, FinancialToolbox, VibratoMonteCarlo, HyperDualNumbers, FinancialFFT, FastGaussQuadrature
#S0 = 100.0;
S0 = taylor_expand(identity, 100.0, order = 5);
K = 105.0;
r = 0.02;
T = 1.0;
d = 0.0;
D = 90.0;

Nsim = 10_000;
Nstep = 30;
sigma = 0.2;
# mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
mc = MonteCarloConfiguration(Nsim, Nstep);

rfCurve = ZeroRate(r);

EUData = EuropeanOption(T, K)
bs = BlackScholesProcess(sigma, Underlying(S0, d));
A = 600.0;
N = 20000;

method = LewisMethod(A, N);

@show result1 = vibrato(bs, rfCurve, mc, EUData, VibratoMonteCarlo.VibratoMonteCarloGaussHermite());
@show result2 = VibratoMonteCarlo.pure_lrm(bs, rfCurve, mc, EUData, VibratoMonteCarlo.VibratoMonteCarloGaussHermite());
@show result3 = blsprice(hyper(S0[0], 1.0, 1.0, 0.0), K, r, T, sigma, d);
@show result4 = [derivative(pricer(bs, rfCurve, method, EUData), i)[0] for i = 0:5];
@show result5 = pricer(bs, rfCurve, mc, EUData);
