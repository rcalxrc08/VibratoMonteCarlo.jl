using TaylorSeries,DualNumbers, FinancialMonteCarlo, FinancialToolbox,VibratoMonteCarlo,HyperDualNumbers,FinancialFFT,FastGaussQuadrature,Test
S0 = 100.0;
# S0 = taylor_expand(identity,S0,order=5);
# S0_d = hyper(S0,1.0,1.0,0.0);
K = 105.0;
r = 0.02;
r_t = taylor_expand(identity,r,order=5);
r_d = hyper(r,1.0,1.0,0.0);
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
rfCurve = ZeroRate(r);

EUData = BinaryEuropeanOption(T, K)
Model = KouProcess(sigma, lam, p, lamp, lamm, Underlying(S0, d));
A = 600.0;
N = 20000;

method = LewisMethod(A, N);

result1 = pure_lrm(Model, ZeroRate(r_t), mc, EUData, VibratoMonteCarlo.VibratoMonteCarloAnalytic(20000,-9.0,19.0));
result2 = pure_lrm(Model, ZeroRate(r_d), mc, EUData, VibratoMonteCarlo.VibratoMonteCarloAnalytic(20000,-9.0,19.0));
right_result = pricer(Model, ZeroRate(r_t), method, EUData)
result4 = pricer(Model, ZeroRate(r_d), method, EUData)
toll=1e-3
right_val=right_result[0];
right_der1=right_result[1];
right_der2=right_result[2]*2;
function check_ders_and_val(val,der1,der2,rv,rd1,rd2)
@test (abs(val-rv)/abs(rv)<toll)
@test (abs(der1-rd1)/abs(rd1)<toll)
@test (abs(der2-rd2)/abs(rd2)<toll)
return;
end
check_ders_and_val(result1[0],result1[1],result1[2]*2,right_val,right_der1,right_der2)
check_ders_and_val(result2.value,result2.epsilon1,result2.epsilon12,right_val,right_der1,right_der2)
check_ders_and_val(result4.value,result4.epsilon1,result4.epsilon12,right_val,right_der1,right_der2)

