using FinancialMonteCarlo, Random;

function generate_last_jump_time_and_value(row_X_i, λ, zero_time, zero_drift, T, mcProcess, dt, mcBaseData)
    t_i = randexp(mcBaseData.parallelMode.rng) / λ
    t_i_1 = t_i / T
    T_jump = zero_time
    x_jump = zero_drift
    last_idx = 0
    jumps_in_step = zero_drift
    while t_i_1 < 1
        jump_size = FinancialMonteCarlo.compute_jump_size(mcProcess, mcBaseData)
        jump_idx = ceil(UInt32, t_i / dt) + true
        @views @. row_X_i[jump_idx:end] += jump_size #add jump component
        jumps_in_step = jump_idx == last_idx ? jumps_in_step + jump_size : jump_size
        last_idx = jump_idx
        # Condition on the grid point before the step containing the last jump, plus the jumps of that step:
        # the diffusion left to T is then exactly Gaussian and lasts at least one step.
        # (The value at jump_idx already contains diffusion beyond t_i, using it with t_i counts that diffusion twice.)
        T_jump = (jump_idx - 2) * dt
        x_jump = row_X_i[jump_idx-1] + jumps_in_step
        t_i += randexp(mcBaseData.parallelMode.rng) / λ
        t_i_1 = t_i / T
    end
    return T_jump, x_jump
end

function simulate_log_with_jump_times(mcProcess::finiteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess <: FinancialMonteCarlo.FiniteActivityProcess, numb <: Number}
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
    S_jumps = Array{typeof(zero_drift)}(undef, Nsim)
    zero_time = 0 * λ
    @inbounds for i = 1:Nsim
        @views row_X_i = X[i, :]
        T_jump, x_jump = generate_last_jump_time_and_value(row_X_i, λ, zero_time, zero_drift, T, mcProcess, dt, mcBaseData)
        @views T_jumps[i] = T_jump
        @views S_jumps[i] = x_jump
    end
    log_S0 = log(S0)
    X_f = @. log_S0 + S_jumps
    return T_jumps, X_f
end
function simulate_with_jump_times(mcProcess::finiteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess <: FinancialMonteCarlo.FiniteActivityProcess, numb <: Number}
    T_jumps, X_f = simulate_log_with_jump_times(mcProcess, rfCurve, mcBaseData, T)
    S_jumps = @. exp(X_f)
    return T_jumps, S_jumps
end