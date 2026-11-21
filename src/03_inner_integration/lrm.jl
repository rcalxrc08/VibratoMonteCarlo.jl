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

# One strike of a smile. Its payoff is FinancialMonteCarlo's payout_strike, the formula a smile and its single options
# share: a smile is integrated strike by strike, without a vector of options per path nor the payoffs of the whole smile
# at every node (#42, #74). Isbits, so it also enters a GPU broadcast.
struct SmileStrike{O, K_}
    K::K_
    isCall::Bool
end
# a scalar in broadcasts, as FinancialMonteCarlo's options: the device kernels broadcast over the paths, not the strike
Base.broadcastable(strike::SmileStrike) = Ref(strike)
# f of every strike of the smile, in a vector: the only allocation of a smile per path
map_strikes(f::F, smile::FinancialMonteCarlo.AbstractSmile{O}) where {F, O} = map((K, isCall) -> f(SmileStrike{O, typeof(K)}(K, isCall)), smile.K, smile.isCall)
payout_f(z, strike::SmileStrike{O}) where {O} = v_value(FinancialMonteCarlo.payout_strike(exp(z), strike.K, strike.isCall, O))
# The density integrates to one for any parameter, so the score terms in v_mod have zero mean and a constant baseline
# (the payoff at the centre of the density, with no derivative) can be subtracted from the payoff and added back:
# the value is unchanged, the derivatives only see payoff - baseline, which removes the baseline times the integration
# error of the score terms (the noise or the bias of the inner integration).
integrand_lrm_base(log_s, log_density_val, eu_opt, baseline) = @muladd (payout_f(log_s, eu_opt) - baseline) * v_mod(log_density_val) + baseline

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
# The samples of a path are drawn once, then every option of a smile is integrated on them, one option at a time: the
# integrand of a whole smile builds its options and two vectors at every sample (12 allocations per sample, #74).
function lrm_interface_vec!(mu, sigma, received_args, mcBaseData, ::VibratoMonteCarloStandard)
    Z, integrand = inner_samples!(received_args.Z, mcBaseData)
    return inner_mean(integrand, Z, mu, sigma, received_args.option)
end

# The inner samples of one path, in the buffer, and the integrand that uses them: antithetic pairs with AntitheticMC and
# with SobolMC (each point paired with its reflection).
inner_samples!(Z, mcBaseData::FinancialMonteCarlo.MonteCarloConfiguration) = (randn!(FinancialMonteCarlo.host_rng(mcBaseData), Z), integrand_lrm_vec)
inner_samples!(Z, mcBaseData::FinancialMonteCarlo.MonteCarloConfiguration{<:FinancialMonteCarlo.AntitheticMC}) = (randn!(FinancialMonteCarlo.host_rng(mcBaseData), Z), integrand_lrm_vec_anti)

# A mean of the integrand at the samples, without an array of its values (#74); one per option for a smile.
inner_mean(integrand::F, Z, mu, sigma, option) where {F} = mean(Base.Fix2(integrand, IntegrandArguments(mu, sigma, option)), Z)
inner_mean(integrand::F, Z, mu, sigma, smile::FinancialMonteCarlo.AbstractSmile) where {F} = map_strikes(strike -> inner_mean(integrand, Z, mu, sigma, strike), smile)

# The trapezoid rule on a uniform grid of n nodes: weights 1, 2, ..., 2, 1, to be multiplied by step / 2. A PatternVector:
# isbits, O(1) memory, on any device and under every AD backend (#69).
trapezoid_weights(n::Integer) = PatternVector(n, PaddedFillPattern(Int8(1), Int8(2), Int8(1)))

# The trapezoid rule on the uniform grid, against the standard normal density (#69), with the weights of init_lrm_vec,
# which include the density: a sum without an array of the values (#74).
function lrm_interface_vec!(mu, sigma, received_args, mcBaseData, ::VibratoMonteCarloAnalytic)
    z, w = received_args.Z
    args = IntegrandArguments(mu, sigma, received_args.option)
    return sum(i -> w[i] * integrand_lrm_vec(z[i], args), eachindex(z, w))
end

function lrm_interface_vec!(mu, sigma, received_args, mcBaseData, ::VibratoMonteCarloGaussHermite)
    Z = received_args.Z
    eu_opt = received_args.option
    x, w = Z
    args = IntegrandArguments(mu, sigma, eu_opt)
    return sum(((xj, wj),) -> integrand_lrm_vec(xj * sqrt2, args) * wj, zip(x, w)) / sqrtπ
end

"""
    singular_points(option)

Log-prices where the payoff of `option` is not smooth (the kink of a call or a put, the jump of a binary), sorted.
VibratoMonteCarloGaussLegendre splits the integral there. Payoffs without a method are integrated as smooth.
"""
singular_points(opt::Union{EuropeanOption, BinaryEuropeanOption, FinancialMonteCarlo.AmericanOption, FinancialMonteCarlo.BinaryAmericanOption}) = (v_value(log(opt.K)),)
singular_points(strike::SmileStrike) = (v_value(log(strike.K)),)
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
    # in the type of sigma, so that Float32 paths stay Float32 (#45)
    z_min = -convert(typeof(v_value_sigma), GAUSS_LEGENDRE_Z)
    z_max = v_value_sigma - z_min
    # The singular points in standard deviations (from values only, so they carry no derivative), still sorted.
    # Outside the bounds the payoff is smooth on the whole range and the piece is skipped.
    cuts = map(p -> clamp((p - v_value_mu) / v_value_sigma, z_min, z_max), singular_points(eu_opt))
    bounds = segment_bounds(z_min, cuts, z_max)
    # A piece of zero length (a cut clamped to a bound) contributes zero, its half-width multiplying every term: no filter,
    # whose sum cannot be inferred. A tuple for an option (its length is known), a vector for a smile.
    return sum(map((za, zb) -> gauss_legendre_segment(za, zb, x, w, args), segment_starts(bounds), segment_ends(bounds)))
end

using Sobol, StatsFuns
# Inner samples from a one dimensional Sobol sequence (cld(Nsim_lrm, 2) uniforms, each point paired with its reflection),
# randomly shifted modulo 1 for every outer path (Cranley-Patterson rotation): each path gets an unbiased estimate
# and their errors are independent.
function init_lrm_vec(x::VibratoMonteCarloStandard, ::FinancialMonteCarlo.SobolConfig, paths)
    seq = SobolSeq(1)
    U = [rng_type(paths)(only(next!(seq))) for _ = 1:cld(x.Nsim_lrm, 2)]
    return (U, similar(U))
end

function inner_samples!((U, Z), mcBaseData::FinancialMonteCarlo.SobolConfig)
    shift = rand(FinancialMonteCarlo.host_rng(mcBaseData))
    @. Z = norminvcdf(mod(U + shift, 1))
    return Z, integrand_lrm_vec_anti
end

# Inner expectation of every outer path, one lrm_interface_vec! per path.
lrm_paths(mu, sigma, args, mcBaseData, vb_mc) = @. lrm_interface_vec!(mu, sigma, args, mcBaseData, vb_mc)

# On the GPU a path cannot run the loops of lrm_interface_vec!, and the configuration (with its CPU rng) cannot enter a
# kernel: the integrand is evaluated on a matrix of paths × nodes in one broadcast, then reduced along the nodes. The
# nodes go to the device with one copy. A smile gives one column per option, on the same inner samples.
lrm_paths(mu::AbstractVector, sigma, args, mcBaseData::CudaConfig, vb_mc) = device_paths(mu, sigma, device_samples(mu, args.Z, mcBaseData, vb_mc), args.option, vb_mc)
device_paths(mu, sigma, samples, option, vb_mc) = device_option_paths(mu, sigma, samples, option, vb_mc)
device_paths(mu, sigma, samples, smile::FinancialMonteCarlo.AbstractSmile, vb_mc) = reduce(hcat, map_strikes(strike -> device_option_paths(mu, sigma, samples, strike, vb_mc), smile))

node_row(mu, x) = permutedims(copyto!(similar(mu, rng_type(mu), length(x)), collect(rng_type(mu), x)))
device_samples(mu, Z, _, ::VibratoMonteCarloGaussLegendre) = (node_row(mu, Z[1]), node_row(mu, Z[2]))
device_samples(mu, (z, w), _, ::VibratoMonteCarloAnalytic) = (node_row(mu, z), node_row(mu, w))
device_samples(mu, Z, _, ::VibratoMonteCarloGaussHermite) = (node_row(mu, Z[1] .* sqrt2), node_row(mu, Z[2]))
# Fresh inner samples for every path, as on the CPU, from the CURAND generator seeded by set_seed!, and antithetic pairs
# with AntitheticMC.
device_samples(mu, Z, mcBaseData, vb_mc::VibratoMonteCarloStandard) = (randn!(similar(mu, rng_type(mu), length(mu), vb_mc.Nsim_lrm)), mcBaseData.sampling isa FinancialMonteCarlo.AntitheticMC)

integrand_at(z, mu, sigma, option) = integrand_lrm_vec(z, IntegrandArguments(mu, sigma, option))
integrand_at_anti(z, mu, sigma, option) = integrand_lrm_vec_anti(z, IntegrandArguments(mu, sigma, option))
weighted_integrand_at(z, mu, sigma, option) = exp(log_density_default(z)) * integrand_at(z, mu, sigma, option)
sum_nodes(f) = vec(sum(f; dims = 2))

function device_option_paths(mu, sigma, (x, w), option, ::VibratoMonteCarloGaussLegendre)
    # The same pieces as lrm_interface_vec!: an empty piece has zero width, so it adds zero instead of being skipped.
    v_mu = v_value.(mu)
    v_sigma = v_value.(sigma)
    z_min = -convert(eltype(v_sigma), GAUSS_LEGENDRE_Z)
    z_max = @. v_sigma - z_min
    cuts = map(p -> (@. clamp((p - v_mu) / v_sigma, z_min, z_max)), singular_points(option))
    # The bounds mix scalars (z_min, and z_max for a scalar sigma) and columns over the paths (the cuts): pieces taken by
    # map over the front and the tail of the tuple, as on the host, are inferred; bounds[i] with i known at run time is not.
    bounds = segment_bounds(z_min, cuts, z_max)
    return sum(map((za, zb) -> sum_nodes(@muladd @. w * (zb - za) / 2 * weighted_integrand_at((za + zb) / 2 + (zb - za) / 2 * x, mu, sigma, option)), segment_starts(bounds), segment_ends(bounds)))
end
device_option_paths(mu, sigma, (z, w), option, ::VibratoMonteCarloAnalytic) = sum_nodes(@. integrand_at(z, mu, sigma, option) * w)
device_option_paths(mu, sigma, (x, w), option, ::VibratoMonteCarloGaussHermite) = sum_nodes(@. integrand_at(x, mu, sigma, option) * w) / sqrtπ
function device_option_paths(mu, sigma, (Z, antithetic), option, ::VibratoMonteCarloStandard)
    f = antithetic ? (@. integrand_at_anti(Z, mu, sigma, option)) : (@. integrand_at(Z, mu, sigma, option))
    return vec(mean(f; dims = 2))
end
