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
# The density is computed in the mode of the configuration: on the GPU in CudaMode, where it stays.
fft_step_law(mcProcess, rfCurve, T, mode = FinancialMonteCarlo.SerialMode()) = fft_step_law(FinancialFFT.density(mcProcess, rfCurve, FinancialMonteCarlo.with_mode(FFTDensityInverter(), mode), T))
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
    density_val = floor_density.(dens, eps_value(first(dens)))
    # Ref: neither the option nor the baseline (a vector for a smile) are broadcast over
    option = Ref(eu_opt)
    baseline = Ref(payout_f(centre, eu_opt))
    f_vals = @. v_value(density_val) * integrand_lrm_base(x, log(density_val), option, baseline)
    # trapezoid: weights 1, 2, ..., 2, 1 times dx / 2
    weights = PatternVector(mc.Npoints, PaddedFillPattern(Int8(1), Int8(2), Int8(1)))
    return law.std * step(z) * sum(weights .* f_vals) / 2
end

# The inner integrals of the paths S_end: one lrm_interface_aug! per path on the host, a vector (of vectors for a smile).
aug_paths(law, option, vb_mc, S_end, mcBaseData) = map(St -> lrm_interface_aug!(law, option, vb_mc, St), S_end)

# On the GPU the density and the paths are on the device: the same trapezoid on a matrix of paths × nodes, evaluated with
# FinancialFFT.density_values (broadcasts and gathers) and reduced along the nodes, as lrm_paths does for the Gaussian
# laws. A smile gives one column per option. density_values builds about ten arrays of the size of the matrix, so the
# paths go in blocks of AUG_BLOCK_NODES nodes (32 MB per array with hyper duals), which keeps the memory bounded.
const AUG_BLOCK_NODES = 2^20
function aug_paths(law::FFTStepLaw, option, mc::VibratoMonteCarloAnalytic, S_end::AbstractVector, ::CudaConfig)
    z = range(mc.x_min, length = mc.Npoints, stop = mc.x_max)
    # trapezoid: weights 1, 2, ..., 2, 1 times dx / 2, with dx = std * step(z)
    trapezoid = PatternVector(mc.Npoints, PaddedFillPattern(1.0, 2.0, 1.0)) .* (law.std * step(z) / 2)
    z_row = node_row(S_end, z)
    w_row = node_row(S_end, trapezoid)
    block = max(1, div(AUG_BLOCK_NODES, mc.Npoints))
    return reduce(vcat, [aug_block(law, option, z_row, w_row, view(S_end, paths)) for paths in Iterators.partition(eachindex(S_end), block)])
end
function aug_block(law::FFTStepLaw, option, z_row, w_row, St)
    log_St = log.(St)
    centre = @. v_value(log_St) + law.mean
    x = @. centre + law.std * z_row
    dens = FinancialFFT.density_values(law.density, x .- log_St)
    return aug_columns(x, centre, dens, w_row, option)
end
aug_columns(x, centre, dens, w_row, option) = (option_ref = Ref(option); density_eps = eps_value(zero(eltype(dens))); sum_nodes(@. aug_integrand(x, centre, dens, density_eps, option_ref) * w_row))
aug_columns(x, centre, dens, w_row, smile::FinancialMonteCarlo.AbstractSmile) = reduce(hcat, [aug_columns(x, centre, dens, w_row, opt) for opt in FinancialMonteCarlo.smile_to_options(smile)])
# The density floored at eps (a float), in the type of the density: max(::Hyper, ::Float64) returns either argument
# unconverted, so it is not inferred (#52), and a broadcast on the GPU needs a concrete type.
floor_density(dens, density_eps) = ifelse(v_value(dens) < density_eps, zero(dens) + density_eps, dens)
# The integrand of lrm_interface_aug! at one node: the density floored at eps, the baseline at the centre.
function aug_integrand(x, centre, dens, density_eps, option)
    density_val = floor_density(dens, density_eps)
    return v_value(density_val) * integrand_lrm_base(x, log(density_val), option, payout_f(centre, option))
end
