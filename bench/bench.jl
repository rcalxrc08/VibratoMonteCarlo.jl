using DualNumbers, HyperDualNumbers, FinancialMonteCarlo, Statistics, BenchmarkTools, VibratoMonteCarlo
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

Nsim = 10_000;
Nsim2 = 20;
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
vb_mc = VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0)
@btime vibrato_saltando($Model, $rfCurve, $mc, $EUData, $vb_mc);
@btime pricer($Model, $rfCurve, $mc, $EUData);

bs = BlackScholesProcess(sigma, Underlying(S0, d));
@btime vibrato($bs, $rfCurve, $mc, $EUData, $vb_mc);
@btime pricer($bs, $rfCurve, $mc, $EUData);

eudatas = [EUData, EUData]
@btime vibrato($bs, $rfCurve, $mc, $eudatas, $vb_mc);
@btime pricer($bs, $rfCurve, $mc, $eudatas);
vb_mc2 = VibratoMonteCarloStandard(Nsim2)
@btime vibrato($bs, $rfCurve, $mc, $eudatas, $vb_mc2);