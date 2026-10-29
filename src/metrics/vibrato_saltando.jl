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
    # X is a Brownian motion subordinated to G: X_T = log(S0) + drift_rn * T + θ * G_T + σ * W(G_T).
    # Given G_T the log price is exactly Gaussian, so the LRM step is taken over the whole subordinator path.
    # Conditioning on the last step only is exact as well, but for VG its increment of G is Gamma with shape dt / κ:
    # it is often so small that σ * sqrt(ΔG) is lost to rounding against log(S) and the LRM weights blow up.
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    θ = mcProcess.θ
    S0 = mcProcess.underlying.S0
    T = get_unique_maturity(eu_opt)
    ψ = compute_drift(mcProcess)
    G_T = rand(mcBaseData.parallelMode.rng, subordinator_law(mcProcess, T), mcBaseData.Nsim)
    drift_rn = r - d - ψ
    mu_jump = @. log(S0) + drift_rn * T + θ * G_T
    sigma_jump = @. σ * sqrt(G_T)
    Z = init_lrm_vec(vb_mc, mcBaseData)
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
    !(typeof(λ) <: AbstractFloat) ? @warn("trying to differentiate over λ, this will lead to incorrect results") : nothing
    !(typeof(T) <: AbstractFloat) ? @warn("trying to differentiate over T, this will lead to incorrect results") : nothing
    (T_jump, S_jump) = VibratoMonteCarlo.simulate_with_jump_times(mcProcess, rfCurve, mcBaseData, T)
    drift_rn = λ * FinancialMonteCarlo.compute_drift(mcProcess)
    result = mean(@. analytic_pricer(S_jump, eu_opt, r, T, T_jump, σ, d, drift_rn) * exp(-r * T_jump))
    return result
end
