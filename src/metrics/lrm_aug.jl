using FinancialFFT, FFTW, Interpolations, PatternVectors
#TODO: switch to FinancialFFT
function characteristic_function_to_density(phi, N, A)
    n = 2^N
    vec = 0.0:(n-1)          # Indices
    dx = 2 * A / n           # Step size, for the density
    x = -A .+ vec * dx       # Grid, for the density
    dt = 2 * pi / (n * dx) # Step size, frequency space
    c = -n / 2 * dt          # Evaluate the characteristic function on [c,d]
    # Both phases, exp(im * vec * dt * A) and exp(-c * x * im), are (-1)^k: dt * A = pi, and n is a multiple of 4.
    # Taking the real part makes the sum a trapezoid: the term at c is the mean of the ones at c and -c, which are conjugate.
    one_minus_one = PatternVector(n, EvenOddPattern(Int8(1), Int8(-1)))
    Y = fft(@. phi(c + vec * dt) * one_minus_one)
    density_vals = @. dt / (2 * pi) * FinancialFFT.real_mod(one_minus_one * Y)
    eps_mod = eps_value(first(density_vals))
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
    spline_cub = cubic_spline_interpolation(x, density_val, extrapolation_bc = Throw())
    w = v_value.(density_val) .* step(x)
    m = sum(x .* w)
    s = sqrt(sum(@. (x - m)^2 * w))
    return DensitySpline(spline_cub, m, s)
end

function adapt_spline(spline_cub, St, x_eval)
    result = @. spline_cub(x_eval - log(St))
    eps_mod = eps_value(first(result))
    y_mod = @. max(result, eps_mod)
    return y_mod
end
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
    weights = PatternVector(mc.Npoints, PaddedFillPattern(Int8(1), Int8(2), Int8(1)))
    return dx * sum(weights .* f_vals) / 2
end
