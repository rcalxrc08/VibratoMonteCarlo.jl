struct LRMInterfaceArguments{num_1, num_2}
    Z::num_1
    option::num_2
    function LRMInterfaceArguments(a::num_1, b::num_2) where {num_1, num_2}
        return new{num_1, num_2}(a, b)
    end
end

Base.broadcastable(x::LRMInterfaceArguments) = Ref(x)
function vibrato(mcProcess::BlackScholesProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    S0 = mcProcess.underlying.S0
    σ = mcProcess.σ
    T = get_unique_maturity(eu_opt)
    drift_rn = (r - d) - σ^2 / 2
    dt = T / mcBaseData.Nstep
    step_vibrato = dt
    X = FinancialMonteCarlo.simulate(BrownianMotion(σ, drift_rn), mcBaseData, T - step_vibrato)
    mu_jump = @views @. X[:, end] + drift_rn * step_vibrato + log(S0)
    sigma_jump = σ * sqrt(step_vibrato)
    Z = init_lrm_vec(vb_mc, mcBaseData)
    #We pack the arguments that could be vectors, but that we don't want to broadcast into a utility struct.
    args = LRMInterfaceArguments(Z, eu_opt)
    # @code_warntype mean(@. lrm_interface_vec!(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)
    result = mean(@. lrm_interface_vec!(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)
    return result
end

function conditional(mcProcess::BlackScholesProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    S0 = mcProcess.underlying.S0
    σ = mcProcess.σ
    T = get_unique_maturity(eu_opt)
    drift_rn = (r - d) - σ^2 / 2
    dt = T / mcBaseData.Nstep
    step_vibrato = dt
    X = FinancialMonteCarlo.simulate(BrownianMotion(σ, drift_rn), mcBaseData, T - step_vibrato)
    mu_jump = @views @. X[:, end] + drift_rn * step_vibrato + log(S0)
    sigma_jump = σ * sqrt(step_vibrato)
    Z = init_lrm_vec(vb_mc, mcBaseData)
    #We pack the arguments that could be vectors, but that we don't want to broadcast into a utility struct.
    args = LRMInterfaceArguments(Z, eu_opt)
    # @code_warntype mean(@. lrm_interface_vec!(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)
    result = mean(@. lrm_interface_vec!(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)
    return result
end

function vibrato(mcProcess::FinancialMonteCarlo.ItoProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    T = get_unique_maturity(eu_opt)
    dt = T / mcBaseData.Nstep
    step_vibrato = dt
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mcBaseData, T - step_vibrato)
    drift_rn = r - d - σ^2 / 2
    mu_jump = @views @. log(S[:, end]) + drift_rn * step_vibrato
    sigma_jump = σ * sqrt(step_vibrato)
    Z = init_lrm_vec(vb_mc, mcBaseData)
    # result = mean(lrm_interface!(Z, mu, sigma_jump, eu_opt, mcBaseData, vb_mc) for mu in mu_jump) * exp(-r * T)
    #We pack the arguments that could be vectors, but that we don't want to broadcast into a utility struct.
    args = LRMInterfaceArguments(Z, eu_opt)
    result = mean(@. lrm_interface_vec!(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)
    return result
end

function vibrato(mcProcess::FinancialMonteCarlo.HestonProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    T = get_unique_maturity(eu_opt)
    dt = T / mcBaseData.Nstep
    (vol, S) = simulate_with_last_vol(mcProcess, rfCurve, mcBaseData, T - dt)
    # if vol is too low, the method does not work.
    drift_rn = r - d
    mu_jump = @views @. log(S[:, end]) + (drift_rn - vol / 2) * dt
    sigma_jump = @. sqrt(dt * vol)
    Z = init_lrm_vec(vb_mc, mcBaseData)
    # Z = init_lrm_vec(vb_mc, mu_jump[1] + sigma_jump[1])
    args = LRMInterfaceArguments(Z, eu_opt)
    result = mean(@. lrm_interface_vec!(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)
    return result
end

function vibrato(mcProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, eu_opt, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    T = get_unique_maturity(eu_opt)
    dt = T / mcBaseData.Nstep
    step_vibrato = dt
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mcBaseData, T - step_vibrato)
    S_end = @views S[:, end]
    adjusted_bound = 18.0
    Z = spline_density(mcProcess, step_vibrato, r, 18, adjusted_bound)
    result = mean(lrm_interface_aug!(mcProcess, Z, eu_opt, mcBaseData, vb_mc, St) for St in S_end) * exp(-r * T)
    return result
end
