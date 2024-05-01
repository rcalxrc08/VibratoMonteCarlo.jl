using FinancialMonteCarlo, Random;
function simulate_with_last_vol(mcProcess::HestonProcess, rfCurve::FinancialMonteCarlo.ZeroRate, mcBaseData::FinancialMonteCarlo.SerialMonteCarloConfig, T::numb) where {numb <: Number}
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = FinancialMonteCarlo.dividend(mcProcess)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    @assert T > 0

    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / κ_s

    dt = T / Nstep
    isDualZeroVol = zero(T * σ₀ * κ * θ * λ1 * σ * ρ)
    isDualZero = isDualZeroVol * S0 * zero(T * r * σ₀ * κ * θ * λ1 * σ * ρ)
    X = Matrix{typeof(isDualZero)}(undef, Nsim, Nstep + 1)
    view(X, :, 1) .= isDualZero
    v_m = [σ₀^2 + isDualZeroVol for _ = 1:Nsim]
    isDualZero_eps = isDualZeroVol + eps(typeof(isDualZeroVol))
    e1 = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim)
    #Wrong, is dual zero shouldn't have rho
    e2_rho = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim)
    e2 = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero) + zero(ρ))}(undef, Nsim)
    tmp_cost = sqrt(1 - ρ^2)
    #TODO: acnaoicna
    for j = 1:Nstep
        randn!(mcBaseData.parallelMode.rng, e1)
        randn!(mcBaseData.parallelMode.rng, e2_rho)
        @. e2 = e1 * ρ + e2_rho * tmp_cost
        @views @. X[:, j+1] = X[:, j] + ((r - d) - v_m / 2) * dt + sqrt(v_m) * sqrt(dt) * e1
        @. v_m += κ_s * (θ_s - v_m) * dt + σ * sqrt(v_m) * sqrt(dt) * e2
        #when v_m = 0.0, the derivative becomes NaN
        @. v_m = max(v_m, isDualZero_eps)
    end
    ## Conclude
    @. X = S0 * exp(X)
    return (v_m, X)
end

function simulate_with_last_vol(mcProcess::HestonProcess, rfCurve::FinancialMonteCarlo.ZeroRate, mcBaseData::FinancialMonteCarlo.SerialAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = FinancialMonteCarlo.dividend(mcProcess)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    @assert T > 0

    ####Simulation
    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / κ_s

    dt = T / Nstep
    isDualZeroVol = zero(T * σ₀ * κ * θ * λ1 * σ * ρ)
    isDualZero = isDualZeroVol * S0 * zero(T * r * σ₀ * κ * θ * λ1 * σ * ρ)
    X = Matrix{typeof(isDualZero)}(undef, Nsim, Nstep + 1)
    view(X, :, 1) .= isDualZero
    isDualZero_eps = isDualZeroVol + eps(typeof(isDualZeroVol))
    Nsim_2 = div(Nsim, 2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]
    v_m = Array{typeof(σ₀ + isDualZeroVol)}(undef, Nsim)
    @. v_m = σ₀^2
    v_m_1 = @views v_m[1:Nsim_2]
    v_m_2 = @views v_m[(Nsim_2+1):end]
    e1 = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim_2)
    e2_rho = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim_2)
    e2 = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, Nsim_2)
    for j = 1:Nstep
        randn!(mcBaseData.parallelMode.rng, e1)
        randn!(mcBaseData.parallelMode.rng, e2_rho)
        @. e2 = -(e1 * ρ + e2_rho * sqrt(1 - ρ * ρ))
        @views @. Xp[:, j+1] = Xp[:, j] + ((r - d) - 0.5 * v_m_1) * dt + sqrt(v_m_1) * sqrt(dt) * e1
        @views @. Xm[:, j+1] = Xm[:, j] + ((r - d) - 0.5 * v_m_2) * dt - sqrt(v_m_2) * sqrt(dt) * e1
        @. v_m_1 += κ_s * (θ_s - v_m_1) * dt + σ * sqrt(v_m_1) * sqrt(dt) * e2
        @. v_m_2 += κ_s * (θ_s - v_m_2) * dt - σ * sqrt(v_m_2) * sqrt(dt) * e2
        @. v_m_1 = max(v_m_1, isDualZero_eps)
        @. v_m_2 = max(v_m_2, isDualZero_eps)
    end
    ## Conclude
    @. X = S0 * exp(X)
    return (v_m, X)
end