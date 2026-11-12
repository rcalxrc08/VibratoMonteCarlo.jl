struct LRMInterfaceArguments{num_1, num_2}
    Z::num_1
    option::num_2
    function LRMInterfaceArguments(a::num_1, b::num_2) where {num_1, num_2}
        return new{num_1, num_2}(a, b)
    end
end

Base.broadcastable(x::LRMInterfaceArguments) = Ref(x)

# A smile with an inner rule on fixed nodes: one integral per option, on the same nodes. The integrand of a whole smile
# builds its options and two vectors at every node (a vector per node and path); Gauss-Legendre also splits each option
# at its own strike only. Not VibratoMonteCarloStandard, which draws new samples at every call.
for inner in (:VibratoMonteCarloGaussLegendre, :VibratoMonteCarloAnalytic, :VibratoMonteCarloGaussHermite)
    @eval function lrm_interface_vec!(mu, sigma, received_args::LRMInterfaceArguments{<:Any, <:FinancialMonteCarlo.AbstractSmile}, mcBaseData, vb_mc::$inner)
        return map(opt -> lrm_interface_vec!(mu, sigma, LRMInterfaceArguments(received_args.Z, opt), mcBaseData, vb_mc), FinancialMonteCarlo.smile_to_options(received_args.option))
    end
end

function pricer(mcProcess::BlackScholesProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::Vibrato, eu_opt_::EuropeanLike)
    # the single options of a smile are built on the host, in every mode: its strikes come to the host first
    eu_opt = FinancialMonteCarlo.to_mode(eu_opt_, FinancialMonteCarlo.SerialMode())
    mcBaseData, vb_mc = method.mc, method.vb_mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    S0 = mcProcess.underlying.S0
    σ = mcProcess.σ
    T = get_unique_maturity(eu_opt)
    drift_rn = (r - d) - drift_correction(mcProcess)
    dt = T / mcBaseData.Nstep
    step_vibrato = dt
    X = FinancialMonteCarlo.simulate(BrownianMotion(σ, drift_rn), mcBaseData, T - step_vibrato)
    mu_jump = @views @muladd @. X[:, end] + drift_rn * step_vibrato + log(S0)
    sigma_jump = σ * sqrt(step_vibrato)
    Z = init_lrm_vec(vb_mc, mcBaseData)
    #We pack the arguments that could be vectors, but that we don't want to broadcast into a utility struct.
    args = LRMInterfaceArguments(Z, eu_opt)
    result = path_mean(lrm_paths(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)
    return result
end

function pricer(mcProcess::BlackScholesProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::Conditional, eu_opt::Union{EuropeanOption, BinaryEuropeanOption})
    mcBaseData = method.mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    T = eu_opt.T
    T_jump = T - T / mcBaseData.Nstep
    S_jump = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mcBaseData, T_jump)
    S_jump_T = @views S_jump[:, end]
    drift_rn = 0
    result = mean(@. analytic_pricer(S_jump_T, eu_opt, r, T, T_jump, σ, d, drift_rn)) * exp(-r * T_jump)
    return result
end

function pricer(mcProcess::FinancialMonteCarlo.ItoProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::Vibrato, eu_opt_::EuropeanLike)
    # the single options of a smile are built on the host, in every mode: its strikes come to the host first
    eu_opt = FinancialMonteCarlo.to_mode(eu_opt_, FinancialMonteCarlo.SerialMode())
    mcBaseData, vb_mc = method.mc, method.vb_mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    T = get_unique_maturity(eu_opt)
    dt = T / mcBaseData.Nstep
    step_vibrato = dt
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mcBaseData, T - step_vibrato)
    drift_rn = r - d - drift_correction(mcProcess)
    mu_jump = @views @muladd @. log(S[:, end]) + drift_rn * step_vibrato
    sigma_jump = σ * sqrt(step_vibrato)
    Z = init_lrm_vec(vb_mc, mcBaseData)
    #We pack the arguments that could be vectors, but that we don't want to broadcast into a utility struct.
    args = LRMInterfaceArguments(Z, eu_opt)
    result = path_mean(lrm_paths(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)
    return result
end

function pricer(mcProcess::FinancialMonteCarlo.HestonProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::Vibrato, eu_opt_::EuropeanLike)
    # the single options of a smile are built on the host, in every mode: its strikes come to the host first
    eu_opt = FinancialMonteCarlo.to_mode(eu_opt_, FinancialMonteCarlo.SerialMode())
    mcBaseData, vb_mc = method.mc, method.vb_mc
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
    result = path_mean(lrm_paths(mu_jump, sigma_jump, args, mcBaseData, vb_mc)) * exp(-r * T)
    return result
end

function pricer(mcProcess::FinancialMonteCarlo.BaseProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::Vibrato, eu_opt_::EuropeanLike)
    # the single options of a smile are built on the host, in every mode: its strikes come to the host first
    eu_opt = FinancialMonteCarlo.to_mode(eu_opt_, FinancialMonteCarlo.SerialMode())
    mcBaseData, vb_mc = method.mc, method.vb_mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    T = get_unique_maturity(eu_opt)
    dt = T / mcBaseData.Nstep
    step_vibrato = dt
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mcBaseData, T - step_vibrato)
    S_end = @views S[:, end]
    # on the device of the paths in CudaMode, the density and the integration too
    law = fft_step_law(mcProcess, rfCurve, step_vibrato, mcBaseData.mode)
    result = path_mean(aug_paths(law, eu_opt, vb_mc, S_end, mcBaseData)) * exp(-r * T)
    return result
end
