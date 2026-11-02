using FinancialMonteCarlo, Random, Distributions;

# The subordinator G_T of VG and NIG, for vibrato_saltando and conditional_saltando: X_T given G_T is Gaussian.

# Stirling series, after shifting the argument above 16 with Γ(x + 1) = x Γ(x). Only elementary functions, because
# DualNumbers, HyperDualNumbers and TaylorSeries do not define loggamma.
@muladd function log_gamma(x)
    shift = max(0, ceil(Int, 16 - v_value(x)))
    x_s = x + shift
    inv_x2 = inv(x_s^2)
    series = evalpoly(inv_x2, (1 / 12, -1 / 360, 1 / 1260, -1 / 1680, 1 / 1188, -691 / 360360, 1 / 156)) / x_s
    return (x_s - 1 / 2) * log(x_s) - x_s + log2π / 2 + series - sum(k -> log(x + k), 0:(shift-1); init = zero(x))
end

# Log density of FinancialMonteCarlo.compute_subordinator_law(mcProcess, κ, t) at g. Distributions.logpdf would need
# loggamma (VG) and is not generic in the number type.
@muladd function subordinator_log_density(::VarianceGammaProcess, g, t, κ)
    shape = t / κ
    return (shape - 1) * log(g) - g / κ - log_gamma(shape) - shape * log(κ)
end
@muladd subordinator_log_density(::NormalInverseGaussianProcess, g, t, κ) = log(t) - (log(κ) + log2π + 3 * log(g) + (g - t)^2 / (κ * g)) / 2

# Draws G_T and its weight. The law depends on κ and T, and a draw is not a smooth function of them (rejection samplers),
# so G_T comes from their values and the dependence goes through the density, in the same v_mod weight of the LRM step.
# Without derivatives on κ and T the weight is true, which leaves the integrand unchanged without computing the density.
function simulate_subordinator(mcProcess, mcBaseData, T)
    κ = mcProcess.κ
    G_T = draw_subordinator(FinancialMonteCarlo.compute_subordinator_law(mcProcess, v_value(κ), v_value(T)), mcBaseData)
    weight = (carries_derivatives(κ) || carries_derivatives(T)) ? map(g -> v_mod(subordinator_log_density(mcProcess, g, T, κ)), G_T) : true
    return (G_T, weight)
end

draw_subordinator(law, mcBaseData) = rand(mcBaseData.parallelMode.rng, law, mcBaseData.Nsim)
# G_T is one dimensional: its quantile at the points of a one dimensional Sobol sequence, skipped as in FinancialMonteCarlo.
function draw_subordinator(law, mcBaseData::FinancialMonteCarlo.SerialSobolMonteCarloConfig)
    seq = SobolSeq(1)
    skip(seq, mcBaseData.Nsim)
    return [quantile(law, only(next!(seq))) for _ in 1:mcBaseData.Nsim]
end

# Mean over the paths with the weight of simulate_subordinator. As in integrand_lrm_base the score terms have zero mean,
# so a baseline with no derivative can be taken out: the sample mean cuts the noise of the κ and T derivatives 3-19 times.
# x has to be discounted already: a discount outside the baseline would multiply its derivative in T by the score noise.
subordinator_mean(x, ::Bool) = mean(x)
function subordinator_mean(x::AbstractVector{<:Number}, weight::AbstractVector)
    baseline = v_value(mean(x))
    return mean(map((x_i, w) -> muladd(x_i - baseline, w, baseline), x, weight))
end
function subordinator_mean(x::AbstractVector{<:AbstractVector}, weight::AbstractVector)
    baseline = v_value.(mean(x))
    return mean(map((x_i, w) -> @.(muladd(x_i - baseline, w, baseline)), x, weight))
end

