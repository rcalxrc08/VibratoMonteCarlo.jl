log_density(x, mu, sigma) = -0.5 * (((x - mu) / sigma)^2 + log(2 * π * sigma^2))
"""
    v_value(x::AbstractFloat)

Returns the value of a generic Number, must be implemented for different backends.
"""
v_value(x::AbstractFloat) = x

"""
    v_mod(x)

Returns the modified values of derivatives in order to automatically implement faa di bruno formula.
"""
function v_mod(x)
    val = exp(v_value(x))
    return exp(x) / val
end

payout_f(z, eu_opt) = FinancialMonteCarlo.payout(exp(z), eu_opt)
integrand_lrm_base(log_s, log_density_val, eu_opt) = v_value(payout_f(log_s, eu_opt)) * v_mod(log_density_val)
integrand_lrm(z, mu, sigma, eu_opt) = integrand_lrm_base(mu + sigma * z, log_density(v_value(mu + sigma * z), mu, sigma), eu_opt)
integrand_lrm_anti(z, mu, sigma, eu_opt) = 0.5 * (integrand_lrm(z, mu, sigma, eu_opt) + integrand_lrm(-z, mu, sigma, eu_opt))

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData::FinancialMonteCarlo.SerialMonteCarloConfig, mc::VibratoMonteCarloStandard)
    # randn!(mcBaseData.parallelMode.rng, Z)
    output = mean(integrand_lrm(z, mu, sigma, eu_opt) for z in Z)
    return output
end

function midpoint_definite_integral(x, f)
    dx = (maximum(x) - minimum(x)) / (length(x) - 1)
    return sum(f(x_) for x_ in x) * dx
end

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData, mc::VibratoMonteCarloAnalytic)
    f(z) = integrand_lrm(z, mu, sigma, eu_opt) * exp(log_density(z, 0, 1))
    return midpoint_definite_integral(Z, f)
end

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData::FinancialMonteCarlo.SerialAntitheticMonteCarloConfig, mc::VibratoMonteCarloStandard)
    # randn!(mcBaseData.parallelMode.rng, Z)
    output = mean(integrand_lrm_anti(z, mu, sigma, eu_opt) for z in Z)
    return output
end

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, ::FinancialMonteCarlo.SerialSobolMonteCarloConfig, mc::VibratoMonteCarloStandard)
    seq = sobolseq(length(Z))
    next!(seq, Z)
    output = mean(integrand_lrm(norminvcdf(z), mu, sigma, eu_opt) for z in Z)
    return output
end
