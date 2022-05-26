log_density(x, mu, sigma) = -0.5 * (((x - mu) / sigma)^2 + log(2 * π * sigma^2))
v_value(x::AbstractFloat) = x
# v_mod(x::AbstractFloat) = 1.0
function v_mod(x)
	val=exp(v_value(x))
	return exp(x)/val;
end

payout_f(z, eu_opt, r) = exp(-r * eu_opt.T) * FinancialMonteCarlo.payout(exp(z), eu_opt)
integrand_lrm_base(log_s, log_density_val, eu_opt, r) = v_value(payout_f(log_s, eu_opt, r)) * v_mod(log_density_val - r * eu_opt.T)
integrand_lrm(z, mu, sigma, eu_opt, r) = integrand_lrm_base(mu + sigma * z, log_density(v_value(mu + sigma * z), mu, sigma), eu_opt, r)
integrand_lrm_anti(z, mu, sigma, eu_opt, r) = 0.5 * (integrand_lrm(z, mu, sigma, eu_opt, r) + integrand_lrm(-z, mu, sigma, eu_opt, r))

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData::FinancialMonteCarlo.SerialMonteCarloConfig, r, mc::VibratoMonteCarloStandard)
    # randn!(mcBaseData.parallelMode.rng, Z)
    output = mean(integrand_lrm(z, mu, sigma, eu_opt, r) for z in Z)
    return output
end

function midpoint_definite_integral(x, f)
    dx = (maximum(x) - minimum(x)) / (length(x) - 1)
    return sum(f(x_) for x_ in x) * dx
end

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData, r, mc::VibratoMonteCarloAnalytic)
    f(z) = integrand_lrm(z, mu, sigma, eu_opt, r) * exp(log_density(z, 0, 1))
    return midpoint_definite_integral(Z, f)
end

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData::FinancialMonteCarlo.SerialAntitheticMonteCarloConfig, r, mc::VibratoMonteCarloStandard)
    # randn!(mcBaseData.parallelMode.rng, Z)
    output = mean(integrand_lrm_anti(z, mu, sigma, eu_opt, r) for z in Z)
    return output
end

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, ::FinancialMonteCarlo.SerialSobolMonteCarloConfig, r, mc::VibratoMonteCarloStandard)
    seq = sobolseq(length(Z))
    next!(seq, Z)
    output = mean(integrand_lrm(norminvcdf(z), mu, sigma, eu_opt, r) for z in Z)
    return output
end
