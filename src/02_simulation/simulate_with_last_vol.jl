using FinancialMonteCarlo, Random;
# Heston paths with the whole variance path: returns (V, S), V[:, j] is the variance at the time of S[:, j]. The same code
# runs on the GPU: the arrays come from FinancialMonteCarlo (CuArrays with CudaMode) and the normals from normals!.
function simulate_with_vol(mcProcess::HestonProcess, und::FinancialMonteCarlo.AbstractUnderlying, rfCurve::FinancialMonteCarlo.ZeroRate{<:Number}, mcBaseData::FinancialMonteCarlo.MonteCarloConfiguration, T::numb) where {numb <: Number}
    r = rfCurve.r
    S0 = und.S0
    d = FinancialMonteCarlo.dividend(und)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    check_domain(T > 0, T, "Time to Maturity must be positive")

    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / κ_s

    dt = T / Nstep
    isDualZeroVol = zero(T * σ₀ * κ * θ * λ1 * σ * ρ)
    isDualZero = isDualZeroVol * S0 * zero(T * r * σ₀ * κ * θ * λ1 * σ * ρ)
    X = FinancialMonteCarlo.get_matrix_type(mcBaseData, mcProcess, typeof(isDualZero))
    view(X, :, 1) .= isDualZero
    V = FinancialMonteCarlo.get_matrix_type(mcBaseData, mcProcess, typeof(σ₀^2 + isDualZeroVol))
    view(V, :, 1) .= σ₀^2 + isDualZeroVol
    isDualZero_eps = isDualZeroVol + eps_value(isDualZeroVol)
    e1 = similar(X, typeof(FinancialMonteCarlo.get_rng_type(isDualZero)), Nsim)
    e2_rho = similar(e1)
    drift, sqrt_dt, ρ_orthogonal = (r - d) * dt, sqrt(dt), sqrt(1 - ρ^2)
    for j = 1:Nstep
        normals!(mcBaseData, e1)
        normals!(mcBaseData, e2_rho)
        # FinancialMonteCarlo's Euler step of the log price and of the variance (fused, with its AD rules, #61); the floor
        # eps_value keeps the variance positive (when v = 0 the derivative of sqrt(v) is NaN)
        heston_columns!(X, V, j, e1, e2_rho, drift, dt, sqrt_dt, κ_s, θ_s, σ, ρ, ρ_orthogonal, isDualZero_eps)
    end
    ## Conclude
    @. X = S0 * exp(X)
    return (V, X)
end

function simulate_with_vol(mcProcess::HestonProcess, und::FinancialMonteCarlo.AbstractUnderlying, rfCurve::FinancialMonteCarlo.ZeroRate{<:Number}, mcBaseData::FinancialMonteCarlo.MonteCarloConfiguration{<:FinancialMonteCarlo.AntitheticMC}, T::numb) where {numb <: Number}
    r = rfCurve.r
    S0 = und.S0
    d = FinancialMonteCarlo.dividend(und)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    check_domain(T > 0, T, "Time to Maturity must be positive")

    ####Simulation
    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / κ_s

    dt = T / Nstep
    isDualZeroVol = zero(T * σ₀ * κ * θ * λ1 * σ * ρ)
    isDualZero = isDualZeroVol * S0 * zero(T * r * σ₀ * κ * θ * λ1 * σ * ρ)
    X = FinancialMonteCarlo.get_matrix_type(mcBaseData, mcProcess, typeof(isDualZero))
    view(X, :, 1) .= isDualZero
    isDualZero_eps = isDualZeroVol + eps_value(isDualZeroVol)
    Nsim_2 = div(Nsim, 2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]
    V = FinancialMonteCarlo.get_matrix_type(mcBaseData, mcProcess, typeof(σ₀^2 + isDualZeroVol))
    view(V, :, 1) .= σ₀^2 + isDualZeroVol
    Vp = @views V[1:Nsim_2, :]
    Vm = @views V[(Nsim_2+1):end, :]
    e1 = similar(X, typeof(FinancialMonteCarlo.get_rng_type(isDualZero)), Nsim_2)
    e2_rho = similar(e1)
    drift, sqrt_dt, ρ_orthogonal = (r - d) * dt, sqrt(dt), sqrt(1 - ρ^2)
    for j = 1:Nstep
        normals!(mcBaseData, e1)
        normals!(mcBaseData, e2_rho)
        # the two halves from the same normals with opposite signs: FinancialMonteCarlo's kernels take the sign in the
        # scale a = ±sqrt(dt) of the increments (#61)
        heston_columns!(Xp, Vp, j, e1, e2_rho, drift, dt, sqrt_dt, κ_s, θ_s, σ, ρ, ρ_orthogonal, isDualZero_eps)
        heston_columns!(Xm, Vm, j, e1, e2_rho, drift, dt, -sqrt_dt, κ_s, θ_s, σ, ρ, ρ_orthogonal, isDualZero_eps)
    end
    ## Conclude
    @. X = S0 * exp(X)
    return (V, X)
end

# Quasi Monte Carlo, as the Sobol Heston of FinancialMonteCarlo: one Sobol point of dimension 2 Nstep per path, the first
# Nstep coordinates drive the price and the others, correlated through ρ, the variance. Without this method a Sobol
# configuration would fall into the MonteCarloConfiguration one above, which matches any sampling, and use the rng.
function simulate_with_vol(mcProcess::HestonProcess, und::FinancialMonteCarlo.AbstractUnderlying, rfCurve::FinancialMonteCarlo.ZeroRate{<:Number}, mcBaseData::FinancialMonteCarlo.SobolConfig, T::numb) where {numb <: Number}
    r = rfCurve.r
    S0 = und.S0
    d = FinancialMonteCarlo.dividend(und)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    check_domain(T > 0, T, "Time to Maturity must be positive")

    κ_s = κ + λ1
    θ_s = κ * θ / κ_s

    dt = T / Nstep
    sqrt_dt = sqrt(dt)
    isDualZeroVol = zero(T * σ₀ * κ * θ * λ1 * σ * ρ)
    isDualZero = isDualZeroVol * S0 * zero(T * r * σ₀ * κ * θ * λ1 * σ * ρ)
    X = Matrix{typeof(isDualZero)}(undef, Nsim, Nstep + 1)
    view(X, :, 1) .= isDualZero
    V = Matrix{typeof(σ₀^2 + isDualZeroVol)}(undef, Nsim, Nstep + 1)
    view(V, :, 1) .= σ₀^2 + isDualZeroVol
    isDualZero_eps = isDualZeroVol + eps_value(isDualZeroVol)
    seq = SobolSeq(2 * Nstep)
    skip(seq, Nstep * Nsim)
    z = Array{typeof(FinancialMonteCarlo.get_rng_type(isDualZero))}(undef, 2 * Nstep)
    ρ_orthogonal = sqrt(1 - ρ^2)
    drift = (r - d) * dt
    for i = 1:Nsim
        next!(seq, z)
        @. z = norminvcdf(z)
        for j = 1:Nstep
            v_m = V[i, j]
            X[i, j+1] = FinancialMonteCarlo.heston_log_step(X[i, j], v_m, z[j], drift, dt, sqrt_dt)
            V[i, j+1] = max(FinancialMonteCarlo.heston_variance_step(v_m, z[j], z[Nstep+j], dt, sqrt_dt, κ_s, θ_s, σ, ρ, ρ_orthogonal), isDualZero_eps)
        end
    end
    ## Conclude
    @. X = S0 * exp(X)
    return (V, X)
end

# One Euler step from column j to column j + 1 of the log prices X and of the variances V, with FinancialMonteCarlo's
# scalar kernels (#61): a = ±sqrt(dt) scales the increments (the sign makes the antithetic half), the variance is floored.
function heston_columns!(X, V, j, e1, e2_rho, drift, dt, a, κ, θ, σ, ρ, ρ_orthogonal, v_floor)
    v = @view V[:, j]
    @views @. X[:, j+1] = FinancialMonteCarlo.heston_log_step(X[:, j], v, e1, drift, dt, a)
    @views @. V[:, j+1] = max(FinancialMonteCarlo.heston_variance_step(v, e1, e2_rho, dt, a, κ, θ, σ, ρ, ρ_orthogonal), v_floor)
    return nothing
end

# Heston paths with the variance at T only.
function simulate_with_last_vol(mcProcess::HestonProcess, und::FinancialMonteCarlo.AbstractUnderlying, rfCurve::FinancialMonteCarlo.ZeroRate{<:Number}, mcBaseData, T::numb) where {numb <: Number}
    (V, X) = simulate_with_vol(mcProcess, und, rfCurve, mcBaseData, T)
    return (V[:, end], X)
end
