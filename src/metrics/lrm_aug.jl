using FinancialFFT, PatternVectors

# The law of one step of the log price for the processes without a closed form density: the FFT density of FinancialFFT
# (FFTDensityInverter), with the mean and the standard deviation of the log increment, values only, used to place the
# integration grid. The density does not depend on S_t: the density of log(S_{t+dt}) given S_t is f(x - log(S_t)).
struct FFTStepLaw{D, F}
    density::D
    mean::F
    std::F
end

function fft_step_law(f::FinancialFFT.FFTDensity)
    # The grid of FinancialFFT is a constant for differentiation (Chernoff range), so a and dx carry no derivative.
    x = range(f.a, step = f.dx, length = length(f.values))
    w = v_value.(f.values) .* f.dx
    m = sum(x .* w)
    s = sqrt(sum(@. (x - m)^2 * w))
    return FFTStepLaw(f, m, s)
end
fft_step_law(mcProcess, rfCurve, T) = fft_step_law(FinancialFFT.density(mcProcess, rfCurve, FFTDensityInverter(), T))
# The same law with values only: the regressions of the American options see no derivative.
value_law(f::FinancialFFT.FFTDensity) = fft_step_law(FinancialFFT.FFTDensity(v_value.(f.values), f.a, f.b, f.dx))

# The inner integral of one path. As for Gaussian processes, x_min and x_max are in standard deviations of the log
# increment around its mean, and the grid is built from values only, so it carries no derivative.
function lrm_interface_aug!(law::FFTStepLaw, eu_opt, mc::VibratoMonteCarloAnalytic, St)
    z = range(mc.x_min, length = mc.Npoints, stop = mc.x_max)
    centre = v_value(log(St)) + law.mean
    # a vector, not a range: shifting a range by an AD number makes a range of AD numbers, which TaylorSeries cannot index
    x = collect(centre .+ law.std .* z)
    # FinancialFFT.density_values is vectorized: all the nodes in one call, shifted by log(St) itself, which carries the
    # pathwise derivatives. The density is zero outside its range: the floor keeps the logarithm finite there.
    dens = FinancialFFT.density_values(law.density, x .- log(St))
    density_val = max.(dens, eps_value(first(dens)))
    # Ref: neither the option nor the baseline (a vector for a smile) are broadcast over
    option = Ref(eu_opt)
    baseline = Ref(payout_f(centre, eu_opt))
    f_vals = @. v_value(density_val) * integrand_lrm_base(x, log(density_val), option, baseline)
    # trapezoid: weights 1, 2, ..., 2, 1 times dx / 2
    weights = PatternVector(mc.Npoints, PaddedFillPattern(Int8(1), Int8(2), Int8(1)))
    return law.std * step(z) * sum(weights .* f_vals) / 2
end
