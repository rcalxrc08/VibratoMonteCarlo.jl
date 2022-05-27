using FinancialMonteCarlo, Random, Distributions;

function simulate_with_sub(mcProcess::SubordinatedBrownianMotion, mcBaseData::FinancialMonteCarlo.SerialMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift * T / mcBaseData.Nstep
    sigma = mcProcess.sigma

    @assert T > 0.0

    type_sub = typeof(rand(mcBaseData.parallelMode.rng, mcProcess.subordinator_))
    isDualZero = drift * zero(type_sub)
    X = Matrix{typeof(isDualZero)}(undef, Nsim, Nstep + 1)
    @views X[:, 1] .= isDualZero
    Z = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim)
    dt_s = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim)
    θ = mcProcess.θ

    @inbounds for i = 1:(Nstep-1)
        rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s)
        randn!(mcBaseData.parallelMode.rng, Z)
        @views @. X[:, i+1] = X[:, i] + drift + θ * dt_s + sigma * sqrt(dt_s) * Z
    end
    # rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s)
    return (dt_s, X)
end

function simulate_with_sub(mcProcess::SubordinatedBrownianMotion, mcBaseData::FinancialMonteCarlo.SerialAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift * T / mcBaseData.Nstep
    sigma = mcProcess.sigma

    @assert T > 0.0

    type_sub = typeof(rand(mcBaseData.parallelMode.rng, mcProcess.subordinator_))
    isDualZero = drift * sigma * zero(type_sub) * 0.0
    X = Matrix{typeof(isDualZero)}(undef, Nsim, Nstep + 1)
    @views X[:, 1] .= isDualZero
    Nsim_2 = div(Nsim, 2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]

    Z = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim_2)
    dt_s = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim_2)
    θ = mcProcess.θ
    @inbounds for j = 1:(Nstep-1)
        rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s)
        randn!(mcBaseData.parallelMode.rng, Z)
        sqrt_dt_s = sqrt.(dt_s)
        @views @. Xp[:, j+1] = Xp[:, j] + dt_s * θ + drift + sqrt_dt_s * Z * sigma
        @views @. Xm[:, j+1] = Xm[:, j] + dt_s * θ + drift - sqrt_dt_s * Z * sigma
    end
    # rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s)
    return (cat(dt_s, dt_s, dims = 1), X)
end

function simulate_with_sub(mcProcess::SubordinatedBrownianMotion, mcBaseData::FinancialMonteCarlo.SerialSobolMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift * T / Nstep
    sigma = mcProcess.sigma
    @assert T > 0
    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    isDualZero = drift * zero(type_sub) * 0.0
    X = Matrix{typeof(isDualZero)}(undef, Nsim, Nstep + 1)
    @views X[:, 1] .= isDualZero
    seq = SobolSeq(Nstep - 1)
    skip(seq, Nstep * Nsim)
    vec = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nstep - 1)
    dt_s_final = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim)
    θ = mcProcess.θ
    @inbounds for i = 1:Nsim
        next!(seq, vec)
        @. vec = norminvcdf(vec)
        dt_s = 0
        @inbounds for j = 1:(Nstep-1)
            dt_s = rand(mcBaseData.parallelMode.rng, mcProcess.subordinator_)
            @views X[i, j+1] = X[i, j] + drift + θ * dt_s + sigma * sqrt(dt_s) * vec[j]
        end
        dt_s_final[i] = dt_s
    end
    return (dt_s_final, X)
end

function compute_drift(mcProcess::VarianceGammaProcess)
    σ = mcProcess.σ
    θ = mcProcess.θ
    κ = mcProcess.κ
    ψ = -1 / κ * log(1 - σ^2 * κ / 2 - θ * κ)
    return ψ
end

function subordinator(mcProcess::VarianceGammaProcess, mcBaseData, T)
    Nstep = mcBaseData.Nstep
    κ = mcProcess.κ
    dt = T / Nstep
    return Gamma(dt / κ, κ)
end

function compute_drift(mcProcess::NormalInverseGaussianProcess)
    σ = mcProcess.σ
    θ = mcProcess.θ
    κ = mcProcess.κ
    ψ = (1 - sqrt(1 - (σ^2 + 2 * θ) * κ)) / κ
    return ψ
end

function subordinator(mcProcess::NormalInverseGaussianProcess, mcBaseData, T)
    Nstep = mcBaseData.Nstep
    κ = mcProcess.κ
    dt = T / Nstep
    return InverseGaussian(dt, (dt^2) / κ)
end

function simulate_with_sub(mcProcess::FinancialMonteCarlo.InfiniteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::Number)
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    @assert T > 0
    ψ = compute_drift(mcProcess)
    drift = r - d - ψ
    θ = mcProcess.θ
    rv = subordinator(mcProcess, mcBaseData, T)

    (dt_s, X) = simulate_with_sub(SubordinatedBrownianMotion(σ, drift, θ, rv), mcBaseData, T)

    S = @. S0 * exp(X)

    return (dt_s, S)
end