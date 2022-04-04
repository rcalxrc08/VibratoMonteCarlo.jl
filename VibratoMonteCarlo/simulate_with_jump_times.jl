using FinancialMonteCarlo, Random;

function simulate_with_jump_times(mcProcess::finiteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess<:FinancialMonteCarlo.FiniteActivityProcess,numb<:Number}

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
    @inbounds for i = 1:Nsim
        t_i = randexp(mcBaseData.parallelMode.rng) / λ
        T_jump = 0.0
        x_jump = zero_drift
        while t_i < T
            T_jump = t_i
            jump_size = FinancialMonteCarlo.compute_jump_size(mcProcess, mcBaseData)
            jump_idx::UInt32 = ceil(UInt32, t_i / dt) + 1
            @views @. X[i, jump_idx:end] += jump_size #add jump component
            x_jump = X[i, jump_idx]
            t_i += randexp(mcBaseData.parallelMode.rng) / λ
        end
		if(T_jump>T-dt)
			@views T_jumps[i] = T_jump
			@views S_jumps[i] = S0 * exp(x_jump)
		else
			@views T_jumps[i] = T-dt
			@views S_jumps[i] = S0 * exp(X[i, end-1])
		end
    end

    return (T_jumps, S_jumps)
end