using DualNumbers, HyperDualNumbers, FinancialMonteCarlo, ForwardDiff, FinancialToolbox
@show "KouModel"
#S0 = 100.0;
S0 = ForwardDiff.Dual(100.0, 0.0, 1.0, 0.0, 0.0);
# S0 = dual(100.0, 1.0);
# S0 = hyper(100.0, 1.0, 1.0, 0.0);
K = 105.0;
#r = 0.02;
r = ForwardDiff.Dual(0.02, 1.0, 0.0, 0.0, 0.0);
# r = hyper(0.02, 1.0, 1.0, 0.0);
# T = ForwardDiff.Dual(1.0,0.0,0.0,1.0,0.0);
T = 1.0;
#T = dual(1.0, 2.0);
d = 0.0;
D = 90.0;

Nsim = 100_000;
Nsim2 = 100;
Nstep = 30;
# sigma = hyper(0.2, 1.0, 1.0, 0.0);
# sigma = dual(0.2, 1.0);
sigma = ForwardDiff.Dual(0.2, 0.0, 0.0, 0.0, 1.0);
p = 0.3;
lam = 5.0;
# lam = dual(5.0, 1.0);
lamp = 30.0;
# lamp = dual(30.0, 1.0);
# lamp = hyper(30.0, 1.0, 1.0, 0.0);
lamm = 20.0;
mu1 = 0.03;
sigma1 = 0.02;
# mc = MonteCarloConfiguration(Nsim, Nstep);
mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
# mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SobolMode());

rfCurve = ZeroRate(r);

EUData = EuropeanOption(T, K)
# Model = KouProcess(sigma, lam, p, lamp, lamm, Underlying(S0, d));
Model = MertonProcess(sigma, lam, mu1, sigma1, Underlying(S0, d));
@show result = vibrato_saltando(Model, rfCurve, mc, EUData, VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0));
@show EuPrice = pricer(Model, rfCurve, mc, EUData);
function blsbin(S0, K, r, T, σ, d)
    d1 = (log(S0 / K) + (r - d + σ * σ * 0.5) * T) / (σ * sqrt(T))
    d2 = d1 - σ * sqrt(T)
    return exp(-r * T) * FinancialToolbox.normcdf(d2)
end

bs = BlackScholesProcess(sigma, Underlying(S0, d));
@show result = vibrato(bs, rfCurve, mc, EUData, VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0));
@show EuPrice = pricer(bs, rfCurve, mc, EUData);
@show collect(ForwardDiff.partials(blsbin(S0, K, r, T, sigma, d)));


# sigma_zero = 0.2;
# kappa = 0.2;
# theta = 0.2;
# lambda = 0.0;
# rho_ = 0.2;
# hest = HestonProcess(sigma, sigma_zero, lambda, kappa, rho_, theta, Underlying(S0, d));
# @show result = vibrato(hest, rfCurve, mc, EUData, VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0));
# @show EuPrice = pricer(hest, rfCurve, mc, EUData);


# vg = VarianceGammaProcess(sigma, 0.01, 0.01, Underlying(S0, d));
# @show result = vibrato_saltando(vg, rfCurve, mc, EUData, VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0));
# @show EuPrice = pricer(vg, rfCurve, mc, EUData);

# nothing