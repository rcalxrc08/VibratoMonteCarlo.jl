function pricer(mcProcess::FinancialMonteCarlo.FiniteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::VibratoSaltando, eu_opt_::EuropeanLike)
    # the single options of a smile are built on the host, in every mode: its strikes come to the host first
    eu_opt = FinancialMonteCarlo.to_mode(eu_opt_, FinancialMonteCarlo.SerialMode())
    mcBaseData, vb_mc = method.mc, method.vb_mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    λ = mcProcess.λ
    T = get_unique_maturity(eu_opt)
    !(typeof(λ) <: AbstractFloat) ? @warn("trying to differentiate over λ, this will lead to incorrect results") : nothing
    !(typeof(T) <: AbstractFloat) ? @warn("trying to differentiate over T, this will lead to incorrect results") : nothing
    (T_jump, X_jump) = simulate_log_with_jump_times(mcProcess, rfCurve, mcBaseData, T)
    drift_rn = (r - d) - drift_correction(mcProcess)
    Z = init_lrm_vec(vb_mc, mcBaseData)
    args = LRMInterfaceArguments(Z, eu_opt)
    mu_jump = @muladd @. X_jump + drift_rn * (T - T_jump)
    sigma_jump = @. σ * sqrt(T - T_jump)
    result = path_mean(lrm_paths(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)

    return result
end

function pricer(mcProcess::FinancialMonteCarlo.InfiniteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::VibratoSaltando, eu_opt_::EuropeanLike)
    # the single options of a smile are built on the host, in every mode: its strikes come to the host first
    eu_opt = FinancialMonteCarlo.to_mode(eu_opt_, FinancialMonteCarlo.SerialMode())
    # X is a Brownian motion subordinated to G: X_T = log(S0) + drift_rn * T + θ * G_T + σ * W(G_T).
    # Given G_T the log price is exactly Gaussian, so the LRM step is taken over the whole subordinator path.
    # Conditioning on the last step only is exact as well, but for VG its increment of G is Gamma with shape dt / κ:
    # it is often so small that σ * sqrt(ΔG) is lost to rounding against log(S) and the LRM weights blow up.
    mcBaseData, vb_mc = method.mc, method.vb_mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    θ = mcProcess.θ
    S0 = mcProcess.underlying.S0
    T = get_unique_maturity(eu_opt)
    (G_T, weight) = simulate_subordinator(mcProcess, mcBaseData, T)
    drift_rn = r - d - drift_correction(mcProcess)
    mu_jump = @muladd @. log(S0) + drift_rn * T + θ * G_T
    sigma_jump = @. σ * sqrt(G_T)
    Z = init_lrm_vec(vb_mc, mcBaseData)
    args = LRMInterfaceArguments(Z, eu_opt)
    discount = exp(-r * T)
    result = subordinator_mean(lrm_paths(mu_jump, sigma_jump, args, mcBaseData, vb_mc) .* discount, weight)
    return result
end

function pricer(mcProcess::FinancialMonteCarlo.FiniteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::ConditionalSaltando, eu_opt::Union{EuropeanOption, BinaryEuropeanOption})
    mcBaseData = method.mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    λ = mcProcess.λ
    T = eu_opt.T
    !(typeof(λ) <: AbstractFloat) ? @warn("trying to differentiate over λ, this will lead to incorrect results") : nothing
    !(typeof(T) <: AbstractFloat) ? @warn("trying to differentiate over T, this will lead to incorrect results") : nothing
    (T_jump, S_jump) = VibratoMonteCarlo.simulate_with_jump_times(mcProcess, rfCurve, mcBaseData, T)
    # The jump part of the drift correction, paid as a dividend yield by the diffusion after the last jump.
    drift_rn = drift_correction(mcProcess) - σ^2 / 2
    result = mean(@. analytic_pricer(S_jump, eu_opt, r, T, T_jump, σ, d, drift_rn) * exp(-r * T_jump))
    return result
end

function pricer(mcProcess::FinancialMonteCarlo.InfiniteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::ConditionalSaltando, eu_opt::Union{EuropeanOption, BinaryEuropeanOption})
    # Given G_T the price is Black-Scholes with total variance σ² * G_T and mean of log(S_T) as in VibratoSaltando,
    # on the same subordinator draws: σ * sqrt(G_T / T) is the volatility over T, the rest of the drift is a dividend yield.
    mcBaseData = method.mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    θ = mcProcess.θ
    S0 = mcProcess.underlying.S0
    T = eu_opt.T
    ψ = drift_correction(mcProcess)
    (G_T, weight) = simulate_subordinator(mcProcess, mcBaseData, T)
    sigma_T = @. σ * sqrt(G_T / T)
    drift_rn = @muladd @. ψ - (θ + σ^2 / 2) * G_T / T
    result = subordinator_mean((@. analytic_pricer(S0, eu_opt, r, T, 0, sigma_T, d, drift_rn)), weight)
    return result
end
