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
v_value(x::Integer) = x

# eps of the float behind x: TaylorSeries does not define eps on its types.
eps_value(x) = eps(typeof(v_value(x)))

# Whether derivatives can flow through x. Mooncake differentiates floats, its extension overlays the method below.
carries_derivatives(::Any) = true
carries_derivatives(::Union{AbstractFloat, Integer}) = false

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
payout_f(z, smile::FinancialMonteCarlo.AbstractSmile) = payout_f_base_v_value.(exp(z), FinancialMonteCarlo.smile_to_options(smile))
# The density integrates to one for any parameter, so the score terms in v_mod have zero mean and a constant baseline
# (the payoff at the centre of the density, with no derivative) can be subtracted from the payoff and added back:
# the value is unchanged, the derivatives only see payoff - baseline, which removes the baseline times the integration
# error of the score terms (the noise or the bias of the inner integration).
integrand_lrm_base(log_s, log_density_val, eu_opt, baseline) = @muladd (payout_f(log_s, eu_opt) - baseline) * v_mod(log_density_val) + baseline
function integrand_lrm_base(log_s, log_density_val, smile::FinancialMonteCarlo.AbstractSmile, baseline)
    payoffs = payout_f(log_s, smile)
    weight = v_mod(log_density_val)
    return @. (payoffs - baseline) * weight + baseline
end

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
# A CudaMode configuration does not land here: every caller goes through lrm_paths, whose device method draws with CURAND.
function lrm_interface_vec!(mu, sigma, received_args, mcBaseData::FinancialMonteCarlo.MonteCarloConfiguration, ::VibratoMonteCarloStandard)
    Z = randn!(FinancialMonteCarlo.host_rng(mcBaseData), received_args.Z)
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

"""
    singular_points(option)

Log-prices where the payoff of `option` is not smooth (the kink of a call or a put, the jump of a binary), sorted.
VibratoMonteCarloGaussLegendre splits the integral there. Payoffs without a method are integrated as smooth.
"""
singular_points(opt::Union{EuropeanOption, BinaryEuropeanOption, FinancialMonteCarlo.AmericanOption, FinancialMonteCarlo.BinaryAmericanOption}) = (v_value(log(opt.K)),)
singular_points(smile::FinancialMonteCarlo.AbstractSmile) = sort!(unique([p for opt in FinancialMonteCarlo.smile_to_options(smile) for p in singular_points(opt)]))
singular_points(::Any) = ()

# ∫_{za}^{zb} of the integrand against the standard normal density, with the Gauss-Legendre nodes x and weights w of [-1, 1].
function gauss_legendre_segment(za, zb, x, w, args)
    c = (za + zb) / 2
    h = (zb - za) / 2
    return sum(begin
        z = @muladd c + h * x[j]
        w[j] * h * exp(log_density_default(z)) * integrand_lrm_vec(z, args)
    end for j in eachindex(x))
end

segment_bounds(z_min, cuts::Tuple, z_max) = (z_min, cuts..., z_max)
segment_bounds(z_min, cuts::AbstractVector, z_max) = [z_min; cuts; z_max]
segment_starts(bounds::Tuple) = Base.front(bounds)
segment_ends(bounds::Tuple) = Base.tail(bounds)
segment_starts(bounds::AbstractVector) = @view bounds[1:(end-1)]
segment_ends(bounds::AbstractVector) = @view bounds[2:end]

function lrm_interface_vec!(mu, sigma, received_args, mcBaseData, ::VibratoMonteCarloGaussLegendre)
    x, w = received_args.Z
    eu_opt = received_args.option
    args = IntegrandArguments(mu, sigma, eu_opt)
    v_value_mu = v_value(mu)
    v_value_sigma = v_value(sigma)
    # A payoff linear in S grows like exp(sigma * z), and exp(sigma * z) * φ(z) = exp(sigma^2 / 2) * φ(z - sigma):
    # its mass is centred at sigma, so the upper bound is moved by sigma (with ±GAUSS_LEGENDRE_Z volga was 4e-9 off at sigma = 0.45).
    z_min = -GAUSS_LEGENDRE_Z
    z_max = GAUSS_LEGENDRE_Z + v_value_sigma
    # The singular points in standard deviations (from values only, so they carry no derivative), still sorted.
    # Outside the bounds the payoff is smooth on the whole range and the piece is skipped.
    cuts = map(p -> clamp((p - v_value_mu) / v_value_sigma, z_min, z_max), singular_points(eu_opt))
    bounds = segment_bounds(z_min, cuts, z_max)
    # A piece of zero length (a cut clamped to a bound) contributes zero, its half-width multiplying every term: no filter,
    # whose sum cannot be inferred. A tuple for an option (its length is known), a vector for a smile.
    return sum(map((za, zb) -> gauss_legendre_segment(za, zb, x, w, args), segment_starts(bounds), segment_ends(bounds)))
end

function lrm_interface_vec!(mu, sigma, received_args, mcBaseData::FinancialMonteCarlo.MonteCarloConfiguration{<:FinancialMonteCarlo.AntitheticMC}, ::VibratoMonteCarloStandard)
    Z = randn!(FinancialMonteCarlo.host_rng(mcBaseData), received_args.Z)
    eu_opt = received_args.option
    args = IntegrandArguments(mu, sigma, eu_opt)
    output = mean(integrand_lrm_vec_anti.(Z, args))
    return output
end
using Sobol, StatsFuns
# Inner samples from a one dimensional Sobol sequence (cld(Nsim_lrm, 2) uniforms, each point paired with its reflection),
# randomly shifted modulo 1 for every outer path (Cranley-Patterson rotation): each path gets an unbiased estimate
# and their errors are independent.
function init_lrm_vec(x::VibratoMonteCarloStandard, ::FinancialMonteCarlo.SobolConfig)
    seq = SobolSeq(1)
    U = [only(next!(seq)) for _ = 1:cld(x.Nsim_lrm, 2)]
    return (U, similar(U))
end

function lrm_interface_vec!(mu, sigma, received_args, mcBaseData::FinancialMonteCarlo.SobolConfig, ::VibratoMonteCarloStandard)
    U, Z = received_args.Z
    eu_opt = received_args.option
    shift = rand(FinancialMonteCarlo.host_rng(mcBaseData))
    @. Z = norminvcdf(mod(U + shift, 1))
    args = IntegrandArguments(mu, sigma, eu_opt)
    output = mean(integrand_lrm_vec_anti.(Z, args))
    return output
end

# Inner expectation of every outer path, one lrm_interface_vec! per path.
lrm_paths(mu, sigma, args, mcBaseData, vb_mc) = @. lrm_interface_vec!(mu, sigma, args, mcBaseData, vb_mc)

# On the GPU a path cannot run the loops of lrm_interface_vec!, and the configuration (with its CPU rng) cannot enter a
# kernel: the integrand is evaluated on a matrix of paths × nodes in one broadcast, then reduced along the nodes. The
# nodes go to the device with one copy. A smile gives one column per option, on the same inner samples.
lrm_paths(mu::AbstractVector, sigma, args, mcBaseData::CudaConfig, vb_mc) = device_paths(mu, sigma, device_samples(mu, args.Z, mcBaseData, vb_mc), args.option, vb_mc)
device_paths(mu, sigma, samples, option, vb_mc) = device_option_paths(mu, sigma, samples, option, vb_mc)
device_paths(mu, sigma, samples, smile::FinancialMonteCarlo.AbstractSmile, vb_mc) = reduce(hcat, [device_option_paths(mu, sigma, samples, opt, vb_mc) for opt in FinancialMonteCarlo.smile_to_options(smile)])

node_row(mu, x) = permutedims(copyto!(similar(mu, Float64, length(x)), collect(Float64, x)))
device_samples(mu, Z, _, ::VibratoMonteCarloGaussLegendre) = (node_row(mu, Z[1]), node_row(mu, Z[2]))
device_samples(mu, Z, _, ::VibratoMonteCarloAnalytic) = (node_row(mu, Z), (maximum(Z) - minimum(Z)) / (length(Z) - 1))
device_samples(mu, Z, _, ::VibratoMonteCarloGaussHermite) = (node_row(mu, Z[1] .* sqrt2), node_row(mu, Z[2]))
# Fresh inner samples for every path, as on the CPU, from the CURAND generator seeded by set_seed!, and antithetic pairs
# with AntitheticMC.
device_samples(mu, Z, mcBaseData, vb_mc::VibratoMonteCarloStandard) = (randn!(similar(mu, Float64, length(mu), vb_mc.Nsim_lrm)), mcBaseData.sampling isa FinancialMonteCarlo.AntitheticMC)

integrand_at(z, mu, sigma, option) = integrand_lrm_vec(z, IntegrandArguments(mu, sigma, option))
integrand_at_anti(z, mu, sigma, option) = integrand_lrm_vec_anti(z, IntegrandArguments(mu, sigma, option))
weighted_integrand_at(z, mu, sigma, option) = exp(log_density_default(z)) * integrand_at(z, mu, sigma, option)
sum_nodes(f) = vec(sum(f; dims = 2))

function device_option_paths(mu, sigma, (x, w), option, ::VibratoMonteCarloGaussLegendre)
    # The same pieces as lrm_interface_vec!: an empty piece has zero width, so it adds zero instead of being skipped.
    option_ref = Ref(option)
    v_mu = v_value.(mu)
    v_sigma = v_value.(sigma)
    z_min = -GAUSS_LEGENDRE_Z
    z_max = @. GAUSS_LEGENDRE_Z + v_sigma
    cuts = map(p -> (@. clamp((p - v_mu) / v_sigma, z_min, z_max)), singular_points(option))
    bounds = (z_min, cuts..., z_max)
    return sum(1:(length(bounds)-1)) do i
        za = bounds[i]
        zb = bounds[i+1]
        sum_nodes(@muladd @. w * (zb - za) / 2 * weighted_integrand_at((za + zb) / 2 + (zb - za) / 2 * x, mu, sigma, option_ref))
    end
end
device_option_paths(mu, sigma, (z, dz), option, ::VibratoMonteCarloAnalytic) = (option_ref = Ref(option); sum_nodes(@. weighted_integrand_at(z, mu, sigma, option_ref)) * dz)
device_option_paths(mu, sigma, (x, w), option, ::VibratoMonteCarloGaussHermite) = (option_ref = Ref(option); sum_nodes(@. integrand_at(x, mu, sigma, option_ref) * w) / sqrtπ)
function device_option_paths(mu, sigma, (Z, antithetic), option, ::VibratoMonteCarloStandard)
    option_ref = Ref(option)
    f = antithetic ? (@. integrand_at_anti(Z, mu, sigma, option_ref)) : (@. integrand_at(Z, mu, sigma, option_ref))
    return vec(mean(f; dims = 2))
end
