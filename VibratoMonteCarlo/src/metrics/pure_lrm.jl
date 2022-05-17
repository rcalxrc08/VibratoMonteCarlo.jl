include("lrm_aug.jl")
# function pure_lrm(mcProcess::FinancialMonteCarlo.ItoProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt::FinancialMonteCarlo.EuropeanPayoff, vb_mc::AbstractVibrato)
    # FinancialMonteCarlo.set_seed!(mcBaseData)
    # r = rfCurve.r
    # d = FinancialMonteCarlo.dividend(mcProcess)
    # σ = mcProcess.σ
    # T = eu_opt.T
    # drift_rn = r - d - σ^2 / 2
	# mu=log(mcProcess.underlying.S0)+drift_rn*T
    # Z = init_lrm_vec(vb_mc, mu)
	# sigma_jump=σ*sqrt(T)
    # return lrm_interface!(Z, mu, sigma_jump, eu_opt, mcBaseData, r,vb_mc)
# end




function pure_lrm(mcProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt::FinancialMonteCarlo.EuropeanPayoff, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
	Z=spline_density(mcProcess,eu_opt.T,r,18,15.0);
    return lrm_interface_aug!(mcProcess,Z, eu_opt, mcBaseData, r,vb_mc)
end