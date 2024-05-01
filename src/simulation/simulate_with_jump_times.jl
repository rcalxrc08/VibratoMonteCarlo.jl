using FinancialMonteCarlo, Random;

function generate_last_jump_time_and_value(row_X_i, λ, zero_time, zero_drift, T, mcProcess, dt, mcBaseData)
    t_i = randexp(mcBaseData.parallelMode.rng) / λ
    T_jump = zero_time
    x_jump = zero_drift
    while t_i < T
        T_jump = t_i
        jump_size = FinancialMonteCarlo.compute_jump_size(mcProcess, mcBaseData)
        jump_idx = ceil(UInt32, t_i / dt) + true
        @views @. row_X_i[jump_idx:end] += jump_size #add jump component
        x_jump = row_X_i[jump_idx]
        t_i += randexp(mcBaseData.parallelMode.rng) / λ
    end
    return T_jump, x_jump
end

function simulate_with_jump_times(mcProcess::finiteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess <: FinancialMonteCarlo.FiniteActivityProcess, numb <: Number}
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    λ = mcProcess.λ

    ####Simulation
    ## Simulate
    # r-d-psi(-i)
    drift_rn = (r - d) - σ^2 / 2 - λ * FinancialMonteCarlo.compute_drift(mcProcess)
    zero_drift = drift_rn * 0 * T
    Nsim = mcBaseData.Nsim
    Nsteps = mcBaseData.Nstep
    X = Matrix{typeof(zero_drift)}(undef, Nsim, Nsteps + 1)
    FinancialMonteCarlo.simulate!(X, BrownianMotion(σ, drift_rn), mcBaseData, T)
    S0 = mcProcess.underlying.S0
    dt = T / Nsteps
    T_jumps = Array{typeof(λ)}(undef, Nsim)
    S_jumps = Array{typeof(zero_drift * S0)}(undef, Nsim)
    zero_time = 0 * λ
    @inbounds for i = 1:Nsim
        @views row_X_i = X[i, :]
        T_jump, x_jump = generate_last_jump_time_and_value(row_X_i, λ, zero_time, zero_drift, T, mcProcess, dt, mcBaseData)
        @views T_jumps[i] = T_jump
        @views S_jumps[i] = S0 * exp(x_jump)
    end

    return T_jumps, S_jumps
end