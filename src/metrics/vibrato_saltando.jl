function vibrato_saltando(mcProcess::finiteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt, vb_mc::AbstractVibrato) where {finiteActivityProcess <: FinancialMonteCarlo.FiniteActivityProcess}
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    λ = mcProcess.λ
    T = get_unique_maturity(eu_opt)
    !(typeof(λ) <: AbstractFloat) ? @warn("trying to differentiate over λ, this will lead to incorrect results") : nothing
    !(typeof(T) <: AbstractFloat) ? @warn("trying to differentiate over T, this will lead to incorrect results") : nothing
    (T_jump, X_jump) = simulate_log_with_jump_times(mcProcess, rfCurve, mcBaseData, T)
    drift_rn = (r - d) - σ^2 / 2 - λ * FinancialMonteCarlo.compute_drift(mcProcess)
    Z = init_lrm_vec(vb_mc, mcBaseData)
    args = LRMInterfaceArguments(Z, eu_opt)
    result = mean(@. lrm_interface_vec!(X_jump + drift_rn * (T - T_jump), σ * sqrt(T - T_jump), args, mcBaseData, vb_mc)) * exp(-r * T)

    return result
end

function vibrato_saltando(mcProcess::infiniteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt, vb_mc::AbstractVibrato) where {infiniteActivityProcess <: FinancialMonteCarlo.InfiniteActivityProcess}
    # Not confirmed theoretically
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    T = get_unique_maturity(eu_opt)
    σ = mcProcess.σ
    ψ = compute_drift(mcProcess)
    (dt_jump, S_jump) = simulate_with_sub(mcProcess, rfCurve, mcBaseData, T)
    drift_rn = r - d - ψ
    #As approximation here I don't consider the law of the last time step in the integration.
    #The last time step is considered the same as the previous, this introduces an error that is proportional to dt.
    mu_jump = @. log(S_jump[:, end-1]) + drift_rn * dt_jump
    sigma_jump = @. σ * sqrt(dt_jump)
    # Z = init_lrm_vec(vb_mc, mu_jump[1])
    Z = init_lrm_vec(vb_mc, mcBaseData)
    # result = mean(lrm_interface!(Z, mu, sigma, eu_opt, mcBaseData, vb_mc) for (mu, sigma) in zip(mu_jump, sigma_jump)) * exp(-r * T)
    args = LRMInterfaceArguments(Z, eu_opt)
    result = mean(@. lrm_interface_vec!(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)
    return result
end

function conditional_saltando(mcProcess::finiteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt) where {finiteActivityProcess <: FinancialMonteCarlo.FiniteActivityProcess}
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    λ = mcProcess.λ
    T = eu_opt.T
    (T_jump, S_jump) = VibratoMonteCarlo.simulate_with_jump_times(mcProcess, rfCurve, mcBaseData, T)
    drift_rn = λ * FinancialMonteCarlo.compute_drift(mcProcess)
    result = mean(@. analytic_pricer(S_jump, eu_opt, r, T, T_jump, σ, d, drift_rn) * exp(-r * T_jump))
    return result
end