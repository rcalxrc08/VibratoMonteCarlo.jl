using FinancialFFT, FFTW, Interpolations
#TODO: switch to FinancialFFT
function characteristic_function_to_density(phi, N, A)
    n = 2^N
    vec = collect(0.0:(n-1))            # Indices
    @views vec[1] = 1e-312
    dx = 2 * A / n           # Step size, for the density
    x = -A .+ (0.0:(n-1)) * dx         # Grid, for the density
    dt = 2 * pi / (n * dx) # Step size, frequency space
    c = -n / 2 * dt          # Evaluate the characteristic function on [c,d]
    # t = @. c + vec * dt         # Grid, frequency space
    Y = @. phi(c + vec * dt) * exp(im * vec * dt * A)
    # @. Y *= exp(im * vec * dt * A)
    fft!(Y)
    density_vals = @. dt / (2 * pi) * FinancialFFT.real_mod(exp(-c * x * im) * Y)
    #TODO: eps should be typed
    eps_mod = eps(eltype(density_vals))
    @. density_vals = max(density_vals, eps_mod)
    return (x, density_vals)
end

function density(mcProcess, T, r, N = 18, xmax = 10.0, St = mcProcess.underlying.S0)
    # characteristic_exponent(x) = FinancialFFT.characteristic_exponent_i(i * x, mcProcess)
    # cf(x) = exp(characteristic_exponent(x) * T)
    corr = FinancialFFT.characteristic_exponent_i(true, mcProcess) * T
    real_cf(x) = exp(((r - mcProcess.underlying.d) * T + log(St) - corr) * im * x + FinancialFFT.characteristic_exponent_i(im * x, mcProcess) * T)
    return characteristic_function_to_density(real_cf, N, xmax)
end

function density_unadapted(mcProcess, T, r, N = 18, xmax = 15.0)
    corr = FinancialFFT.characteristic_exponent_i(true, mcProcess) * T
    real_cf(x) = exp(((r - mcProcess.underlying.d) * T - corr) * 1im * x + FinancialFFT.characteristic_exponent_i(im * x, mcProcess) * T)
    return characteristic_function_to_density(real_cf, N, xmax)
end

# Cubic spline of the density of the log increment, with its mean and standard deviation (values only, no derivatives)
# used to place the integration grid. Calling it evaluates the spline.
struct DensitySpline{S, F}
    spline::S
    mean::F
    std::F
end
(d::DensitySpline)(y) = d.spline(y)

spline_density(mcProcess, T, r, N = 18, xmax = 15.0) = density_spline(density_unadapted(mcProcess, T, r, N, xmax)...)

function density_spline(x, density_val)
    spline_cub = CubicSplineInterpolation(x, density_val, extrapolation_bc = Throw())
    w = v_value.(density_val) .* step(x)
    m = sum(x .* w)
    s = sqrt(sum(@. (x - m)^2 * w))
    return DensitySpline(spline_cub, m, s)
end

function adapt_spline(spline_cub, St, x_eval)
    result = @. spline_cub(x_eval - log(St))
    #TODO: eps should be typed
    eps_mod = eps(eltype(result))
    y_mod = @. max(result, eps_mod)
    return y_mod
end
#TODO: use AlternateVectors
function lrm_interface_aug!(mcProcess, spline_dens::DensitySpline, eu_opt, mcBaseData, mc::VibratoMonteCarloAnalytic, St = mcProcess.underlying.S0)
    # As for Gaussian processes, x_min and x_max are in standard deviations of the log increment around its mean.
    # The grid is built from values only, so it carries no derivative.
    z = range(mc.x_min, length = mc.Npoints, stop = mc.x_max)
    centre = v_value(log(St)) + spline_dens.mean
    x = @. centre + spline_dens.std * z
    density_val = adapt_spline(spline_dens, St, x)
    dx = spline_dens.std * step(z)
    baseline = payout_f(centre, eu_opt)
    f_vals = [v_value(dens) * integrand_lrm_base(log_s, log(dens), eu_opt, baseline) for (log_s, dens) in zip(x, density_val)]
    @views result = dx * (f_vals[1] + f_vals[end] + 2 * sum(f_vals[2:(end-1)])) / 2
    return result
end
