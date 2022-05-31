function vibrato(mcProcess::FinancialMonteCarlo.ItoProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt::FinancialMonteCarlo.EuropeanPayoff, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    T = eu_opt.T
    dt = T / mcBaseData.Nstep
    step_vibrato = dt
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mcBaseData, T - step_vibrato)
    drift_rn = r - d - σ^2 / 2
    mu_jump = @views @. log(S[:, end]) + drift_rn * step_vibrato
    sigma_jump = σ * sqrt(step_vibrato)
    Z = init_lrm_vec(vb_mc, mcBaseData)
    result = mean(lrm_interface!(Z, mu, sigma_jump, eu_opt, mcBaseData, vb_mc) for mu in mu_jump)*exp(-r*T)
    return result
end

function vibrato(mcProcess::FinancialMonteCarlo.HestonProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt::FinancialMonteCarlo.EuropeanPayoff, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    T = eu_opt.T
    dt = T / mcBaseData.Nstep
    (vol, S) = simulate_with_last_vol(mcProcess, rfCurve, mcBaseData, T - dt)
    # if vol is too low, the method does not work.
    drift_rn = r - d
    mu_jump = @views @. log(S[:, end]) + (drift_rn - 0.5 * vol) * dt
    sigma_jump = @. sqrt(dt * vol)
    Z = init_lrm_vec(vb_mc, mu_jump[1] + sigma_jump[1])
    result = mean(lrm_interface!(Z, mu, sigma, eu_opt, mcBaseData, vb_mc) for (mu, sigma) in zip(mu_jump, sigma_jump))*exp(-r*T)
    return result
end

function vibrato(mcProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt::FinancialMonteCarlo.EuropeanPayoff, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    T = eu_opt.T
    dt = T / mcBaseData.Nstep
    step_vibrato = dt
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mcBaseData, T - step_vibrato)
    S_end = @views S[:, end]
    Z = spline_density(mcProcess, step_vibrato, r, 18, 20.0)
    result = mean(lrm_interface_aug!(mcProcess, Z, eu_opt, mcBaseData, vb_mc, St) for St in S_end)*exp(-r*T)
    return result
end

function vibrato(mcProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opts, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    T = maximum(eu_opt.T for eu_opt in eu_opts)
    dt = T / mcBaseData.Nstep
    step_vibrato = dt
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mcBaseData, T - step_vibrato)
    S_end = @views S[:, end]
    Z = spline_density(mcProcess, step_vibrato, r, 18, 20.0)
    result = [mean(lrm_interface_aug!.(mcProcess, Z, eu_opt, mcBaseData, vb_mc, St)) for St in S_end]
    return result
end