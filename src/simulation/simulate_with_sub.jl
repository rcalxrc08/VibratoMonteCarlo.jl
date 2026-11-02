using FinancialMonteCarlo, Random, Distributions;

function simulate_with_sub(mcProcess::SubordinatedBrownianMotion, mcBaseData::FinancialMonteCarlo.SerialMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift * T / mcBaseData.Nstep
    sigma = mcProcess.sigma

    @assert T > 0

    type_sub = typeof(rand(mcBaseData.parallelMode.rng, mcProcess.subordinator_))
    isDualZero = drift * zero(type_sub)
    X = Matrix{typeof(isDualZero)}(undef, Nsim, Nstep + 1)
    @views X[:, 1] .= isDualZero
    Z = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim)
    dt_s = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim)
    θ = mcProcess.θ

    @inbounds for i = 1:Nstep
        rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s)
        randn!(mcBaseData.parallelMode.rng, Z)
        @views @. X[:, i+1] = X[:, i] + drift + θ * dt_s + sigma * sqrt(dt_s) * Z
    end
    return (dt_s, X)
end

function simulate_with_sub(mcProcess::SubordinatedBrownianMotion, mcBaseData::FinancialMonteCarlo.SerialAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift * T / mcBaseData.Nstep
    sigma = mcProcess.sigma

    @assert T > 0

    type_sub = typeof(rand(mcBaseData.parallelMode.rng, mcProcess.subordinator_))
    isDualZero = drift * sigma * zero(type_sub) * 0
    X = Matrix{typeof(isDualZero)}(undef, Nsim, Nstep + 1)
    @views X[:, 1] .= isDualZero
    Nsim_2 = div(Nsim, 2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]

    Z = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim_2)
    dt_s = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim_2)
    θ = mcProcess.θ
    @inbounds for j = 1:Nstep
        rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s)
        randn!(mcBaseData.parallelMode.rng, Z)
        sqrt_dt_s = sqrt.(dt_s)
        @views @. Xp[:, j+1] = Xp[:, j] + dt_s * θ + drift + sqrt_dt_s * Z * sigma
        @views @. Xm[:, j+1] = Xm[:, j] + dt_s * θ + drift - sqrt_dt_s * Z * sigma
    end
    return (cat(dt_s, dt_s, dims = 1), X)
end

function simulate_with_sub(mcProcess::SubordinatedBrownianMotion, mcBaseData::FinancialMonteCarlo.SerialSobolMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift * T / Nstep
    sigma = mcProcess.sigma
    @assert T > 0
    #TODO: adjust type
    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    isDualZero = drift * zero(type_sub)
    X = Matrix{typeof(isDualZero)}(undef, Nsim, Nstep + 1)
    @views X[:, 1] .= isDualZero
    seq = SobolSeq(Nstep)
    skip(seq, Nstep * Nsim)
    vec = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nstep)
    dt_s_final = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim)
    θ = mcProcess.θ
    @inbounds for i = 1:Nsim
        next!(seq, vec)
        @. vec = norminvcdf(vec)
        dt_s = 0
        @inbounds for j = 1:Nstep
            dt_s = rand(mcBaseData.parallelMode.rng, mcProcess.subordinator_)
            @views X[i, j+1] = X[i, j] + drift + θ * dt_s + sigma * sqrt(dt_s) * vec[j]
        end
        dt_s_final[i] = dt_s
    end
    return (dt_s_final, X)
end

# Law of the subordinator increment over a time interval t (mean t, variance κ t), from FinancialMonteCarlo.
subordinator(mcProcess, mcBaseData, T) = FinancialMonteCarlo.compute_subordinator_law(mcProcess, T / mcBaseData.Nstep)

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
    G_T = rand(mcBaseData.parallelMode.rng, FinancialMonteCarlo.compute_subordinator_law(mcProcess, v_value(κ), v_value(T)), mcBaseData.Nsim)
    weight = (carries_derivatives(κ) || carries_derivatives(T)) ? map(g -> v_mod(subordinator_log_density(mcProcess, g, T, κ)), G_T) : true
    return (G_T, weight)
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

function simulate_with_sub(mcProcess::FinancialMonteCarlo.InfiniteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::Number)
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    @assert T > 0
    drift = r - d - drift_correction(mcProcess)
    θ = mcProcess.θ
    rv = subordinator(mcProcess, mcBaseData, T)

    (dt_s, X) = simulate_with_sub(SubordinatedBrownianMotion(σ, drift, θ, rv), mcBaseData, T)

    S = @. S0 * exp(X)

    return (dt_s, S)
end