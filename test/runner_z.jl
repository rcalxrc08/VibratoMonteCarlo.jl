S0 = 100.0;
K = 80.0;
r = 0.01;
T = 1.0;
d = 0.01;
D = 90.0;

Nsim = 10_000;
Nsim2 = 20;
Nstep = 30;
sigma = 0.2;
mc = MonteCarloConfiguration(Nsim, Nstep);

@show result = gradient((S0,K,r,T,sigma) -> vibrato_z(BlackScholesProcess(sigma, Underlying(S0, d)), ZeroRate(r), mc, EuropeanOption(T, K), VibratoMonteCarloAnalytic(Nsim2, -5.0, 5.0)),S0,K,r,T,sigma);

nothing