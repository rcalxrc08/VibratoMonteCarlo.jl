function pricer(mcProcess::FinancialMonteCarlo.ItoProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::PureLRM, eu_opt_::EuropeanLike)
    # the single options of a smile are built on the host, in every mode: its strikes come to the host first
    eu_opt = FinancialMonteCarlo.to_mode(eu_opt_, FinancialMonteCarlo.SerialMode())
    mcBaseData, vb_mc = method.mc, method.vb_mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    T = get_unique_maturity(eu_opt)
    drift_rn = r - d - drift_correction(mcProcess)
    mu = @muladd log(mcProcess.underlying.S0) + drift_rn * T
    Z = init_lrm_vec(vb_mc, mcBaseData)
    sigma_jump = σ * sqrt(T)
    args = LRMInterfaceArguments(Z, eu_opt)
    result = lrm_interface_vec!(mu, sigma_jump, args, mcBaseData, vb_mc) * exp(-r * T)
    return result
end

function pricer(mcProcess::FinancialMonteCarlo.BaseProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::PureLRM, eu_opt_::EuropeanLike)
    # the single options of a smile are built on the host, in every mode: its strikes come to the host first
    eu_opt = FinancialMonteCarlo.to_mode(eu_opt_, FinancialMonteCarlo.SerialMode())
    mcBaseData, vb_mc = method.mc, method.vb_mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    T = get_unique_maturity(eu_opt)
    Z = fft_step_law(mcProcess, rfCurve, T)
    return lrm_interface_aug!(Z, eu_opt, vb_mc, mcProcess.underlying.S0) * exp(-r * T)
end