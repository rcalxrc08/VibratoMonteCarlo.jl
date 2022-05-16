payout_f_aug(z, eu_opt, r) = exp(-r * eu_opt.T) * FinancialMonteCarlo.payout(exp(z), eu_opt)
using FinancialFFT,FFTW,Interpolations
function characteristic_function_to_density(phi,N,A) 
  n=2^N;
  vec = collect(0.0:(n-1))            # Indices
  @views vec[1]=1e-312;
  dx = 2*A/n           # Step size, for the density
  x = -A .+ (0.0:(n-1)) * dx         # Grid, for the density
  dt = 2*pi / ( n * dx ) # Step size, frequency space
  c = -n/2 * dt          # Evaluate the characteristic function on [c,d]
  t = c .+ vec * dt         # Grid, frequency space
  phi_t = phi.(t)
  Y = exp.( 1im * vec * dt * A ) .* phi_t
  fft!(Y)
  density_vals = max.(FinancialFFT.real_mod.(dt / (2*pi) * exp.( - 1im * c * x ) .* Y),eps())
  return (x,density_vals);
end

function density(mcProcess,T,r,N=18,xmax=10.0,St=mcProcess.underlying.S0) 
	cf = FinancialFFT.CharactheristicFunction(mcProcess, T)
	corr = FinancialFFT.CharactheristicExponent(-1im, mcProcess, T)
	real_cf(x)=cf(x)*exp(((r-mcProcess.underlying.d)*T+log(St)-corr )* 1im * x)
	return characteristic_function_to_density(real_cf,N,xmax)
end

function density2(mcProcess,T,r,N=18,xmax=15.0) 
	cf = FinancialFFT.CharactheristicFunction(mcProcess, T)
	corr = FinancialFFT.CharactheristicExponent(-1im, mcProcess, T)
	real_cf(x)=cf(x)*exp(((r-mcProcess.underlying.d)*T-corr )* 1im * x)
	return characteristic_function_to_density(real_cf,N,xmax)
end


function adapt_density(x,density_val,St,x_eval)
	spline_cub = CubicSplineInterpolation(x, density_val,extrapolation_bc=Flat())
	y_mod=max.(spline_cub.(x_eval.-log(St)),eps())
	return (x_eval,y_mod)
end

function spline_density(mcProcess,T,r,N=18,xmax=15.0)
	(x,density_val)=density2(mcProcess,T,r,N,xmax);
	spline_cub = CubicSplineInterpolation(x, density_val,extrapolation_bc=Flat())
	return spline_cub
end

function adapt_spline(spline_cub,St,x_eval)
	y_mod=max.(spline_cub.(x_eval.-log(St)),eps())
	return y_mod
end

function lrm_interface_aug!(mcProcess,Z, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData, r,mc::VibratoMonteCarloAnalytic,St=mcProcess.underlying.S0)
	spline_dens=Z;
	x=range(mc.x_min, length=mc.Npoints, stop=mc.x_max)
	density_val=adapt_spline(spline_dens,St,x);
	dx=x.step.hi;
	f_vals=[v_value(d*payout_f_aug(x_v, eu_opt, r)) *  v_mod(log(d) - r * eu_opt.T) for (x_v,d) in zip(x,density_val)];
	@views result=0.5*dx*(f_vals[1]+f_vals[end]+ 2.0*sum(f_vals[2:(end-1)]))
	return result
end

