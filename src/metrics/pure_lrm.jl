function pricer(mcProcess::FinancialMonteCarlo.ItoProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::PureLRM, eu_opt::EuropeanLike)
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

function pricer(mcProcess::FinancialMonteCarlo.BaseProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::PureLRM, eu_opt::EuropeanLike)
    mcBaseData, vb_mc = method.mc, method.vb_mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    T = get_unique_maturity(eu_opt)
    Z = spline_density(mcProcess, T, r, 18, 15.0)
    return lrm_interface_aug!(mcProcess, Z, eu_opt, mcBaseData, vb_mc) * exp(-r * T)
end