function log_density_default(x)
    return (x^2 + log2π) / (-2)
end
function log_density(x_in, mu, sigma)
    x = (x_in - mu) / sigma
    return log_density_default(x) - log(sigma)
end
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
    return exp(x - v_value(x))
end
v_mod(::AbstractFloat) = 1

payout_f_base_v_value(s, eu_opt) = v_value(FinancialMonteCarlo.payout(s, eu_opt))

payout_f(z, eu_opt) = payout_f_base_v_value(exp(z), eu_opt)
function payout_f(z, eu_opts::Array)
    S = exp(z)
    return payout_f_base_v_value.(S, eu_opts)
end
# The density integrates to one for any parameter, so the score terms in v_mod have zero mean and a constant baseline
# (the payoff at the centre of the density, with no derivative) can be subtracted from the payoff and added back:
# the value is unchanged, the derivatives only see payoff - baseline, which removes the baseline times the integration
# error of the score terms (the noise or the bias of the inner integration).
integrand_lrm_base(log_s, log_density_val, eu_opt, baseline) = @muladd (payout_f(log_s, eu_opt) - baseline) * v_mod(log_density_val) + baseline
integrand_lrm_base(log_s, log_density_val, eu_opts::Array, baseline) = @. (payout_f(log_s, eu_opts) - baseline) * v_mod(log_density_val) + baseline

struct IntegrandArguments{num_1, num_2, num_3}
    mu::num_1
    sigma::num_2
    option::num_3
    function IntegrandArguments(a::num_1, b::num_2, c::num_3) where {num_1, num_2, num_3}
        return new{num_1, num_2, num_3}(a, b, c)
    end
end

Base.broadcastable(x::IntegrandArguments) = Ref(x)

function integrand_lrm_vec(z, args)
    mu = args.mu
    sigma = args.sigma
    eu_opt = args.option
    v_value_mu = v_value(mu)
    z_adj = @muladd v_value_mu + v_value(sigma) * z
    return integrand_lrm_base(z_adj, log_density(z_adj, mu, sigma), eu_opt, payout_f(v_value_mu, eu_opt))
end

function integrand_lrm_vec_anti(z, args)
    mu = args.mu
    sigma = args.sigma
    eu_opt = args.option
    v_value_mu = v_value(mu)
    v_value_sigma = v_value(sigma)
    baseline = payout_f(v_value_mu, eu_opt)
    z_scaled = v_value_sigma * z
    z_adj_p = v_value_mu + z_scaled
    res_p = integrand_lrm_base(z_adj_p, log_density(z_adj_p, mu, sigma), eu_opt, baseline)
    z_adj_m = v_value_mu - z_scaled
    res_m = integrand_lrm_base(z_adj_m, log_density(z_adj_m, mu, sigma), eu_opt, baseline)
    return (res_m + res_p) / 2
end

# With VibratoMonteCarloStandard the inner samples are drawn again for every outer path (Z is only a buffer):
# sharing them would turn their error into a bias that does not average out over the outer paths.
function lrm_interface_vec!(mu, sigma, received_args, mcBaseData::FinancialMonteCarlo.SerialMonteCarloConfig, ::VibratoMonteCarloStandard)
    Z = randn!(mcBaseData.parallelMode.rng, received_args.Z)
    eu_opt = received_args.option
    args = IntegrandArguments(mu, sigma, eu_opt)
    output = mean(integrand_lrm_vec.(Z, args))
    return output
end

function midpoint_definite_integral_vec(x, f)
    dx = (maximum(x) - minimum(x)) / (length(x) - 1)
    return sum(f) * dx
end

function lrm_interface_vec!(mu, sigma, received_args, mcBaseData, ::VibratoMonteCarloAnalytic)
    Z = received_args.Z
    eu_opt = received_args.option
    args = IntegrandArguments(mu, sigma, eu_opt)
    f = @. integrand_lrm_vec(Z, args) * exp(log_density_default(Z))
    return midpoint_definite_integral_vec(Z, f)
end

function lrm_interface_vec!(mu, sigma, received_args, mcBaseData, ::VibratoMonteCarloGaussHermite)
    Z = received_args.Z
    eu_opt = received_args.option
    x, w = Z
    args = IntegrandArguments(mu, sigma, eu_opt)
    output = sum(@. integrand_lrm_vec(x * sqrt2, args) * w) / sqrtπ
    return output
end

function lrm_interface_vec!(mu, sigma, received_args, mcBaseData::FinancialMonteCarlo.SerialAntitheticMonteCarloConfig, ::VibratoMonteCarloStandard)
    Z = randn!(mcBaseData.parallelMode.rng, received_args.Z)
    eu_opt = received_args.option
    args = IntegrandArguments(mu, sigma, eu_opt)
    output = mean(integrand_lrm_vec_anti.(Z, args))
    return output
end
using Sobol, StatsFuns
# Inner samples from a one dimensional Sobol sequence (cld(Nsim_lrm, 2) uniforms, each point paired with its reflection),
# randomly shifted modulo 1 for every outer path (Cranley-Patterson rotation): each path gets an unbiased estimate
# and their errors are independent.
function init_lrm_vec(x::VibratoMonteCarloStandard, ::FinancialMonteCarlo.SerialSobolMonteCarloConfig)
    seq = SobolSeq(1)
    U = [only(next!(seq)) for _ = 1:cld(x.Nsim_lrm, 2)]
    return (U, similar(U))
end

function lrm_interface_vec!(mu, sigma, received_args, mcBaseData::FinancialMonteCarlo.SerialSobolMonteCarloConfig, ::VibratoMonteCarloStandard)
    U, Z = received_args.Z
    eu_opt = received_args.option
    shift = rand(mcBaseData.parallelMode.rng)
    @. Z = norminvcdf(mod(U + shift, 1))
    args = IntegrandArguments(mu, sigma, eu_opt)
    output = mean(integrand_lrm_vec_anti.(Z, args))
    return output
end
