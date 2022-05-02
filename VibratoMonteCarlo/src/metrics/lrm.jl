log_density(x, mu, sigma) = -0.5 * (((x - mu) / sigma)^2 + log(2 * π * sigma^2))
v_value(x::AbstractFloat) = x

payout_f(z, mu, sigma, eu_opt, r) = exp(-r * eu_opt.T) * FinancialMonteCarlo.payout(exp(mu + sigma * z), eu_opt)
integrand_lrm(z, mu, sigma, eu_opt, r) = v_value(payout_f(z, mu, sigma, eu_opt, r)) * v_mod(log_density(v_value(mu + sigma * z), mu, sigma) - r * eu_opt.T)
integrand_lrm_anti(z, mu, sigma, eu_opt, r) = 0.5 * (integrand_lrm(z, mu, sigma, eu_opt, r) + integrand_lrm(-z, mu, sigma, eu_opt, r))

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData::FinancialMonteCarlo.SerialMonteCarloConfig, r,mc::VibratoMonteCarloStandard)
    randn!(mcBaseData.parallelMode.rng, Z)
    output = mean(integrand_lrm(z, mu, sigma, eu_opt, r) for z in Z)
    return output
end

function midpoint_definite_integral(x, f)
    dx = (maximum(x) - minimum(x)) / (length(x) - 1)
    return reduce(+, f(x_) * dx for x_ in x)
end

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData, r,mc::VibratoMonteCarloAnalytic)
    f(z) = integrand_lrm(z, mu, sigma, eu_opt, r) * exp(log_density(z, 0, 1))
    return midpoint_definite_integral(Z, f)
end


function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData::FinancialMonteCarlo.SerialAntitheticMonteCarloConfig, r,mc::VibratoMonteCarloStandard)
    randn!(mcBaseData.parallelMode.rng, Z)
    output = mean(integrand_lrm_anti(z, mu, sigma, eu_opt, r) for z in Z)
    return output
end

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, ::FinancialMonteCarlo.SerialSobolMonteCarloConfig, r,mc::VibratoMonteCarloStandard)
    seq = sobolseq(length(Z))
    next!(seq, Z)
    output = mean(integrand_lrm(norminvcdf(z), mu, sigma, eu_opt, r) for z in Z)
    return output
end

