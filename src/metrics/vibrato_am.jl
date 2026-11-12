using LinearAlgebra

# Continuation value regressed on 1, x, x^2 (Longstaff-Schwartz), values only.
continuation_coefficients(x, y) = hcat(one.(x), x, x .^ 2) \ y
continuation_value(α, x) = α[1] + x * (α[2] + x * α[3])

# American options: see the docstring of Vibrato.
# The regressions run on the host on the values of the paths, taken there once; the smoothed payoffs (the vibrato
# expectations over one step, values only for the regressions, with derivatives for the cash flows) are computed for a
# set of paths at a time by step_payoffs, on the device of the paths in CudaMode.
function pricer(mcProcess::FinancialMonteCarlo.BaseProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::Vibrato, am_opt::FinancialMonteCarlo.AmericanPayoff)
    mcBaseData, vb_mc = method.mc, method.vb_mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    T = am_opt.T
    Nstep = mcBaseData.Nstep
    dt = T / Nstep
    S, laws = frozen_laws(mcProcess, rfCurve, mcBaseData, T)
    S_value = v_value.(S)
    S_value_host = to_host(S_value)
    Nsim = size(S, 1)
    args = LRMInterfaceArguments(init_lrm_vec(vb_mc, mcBaseData), am_opt)
    # the smoothed values of the paths `paths` at the column k, on the host
    smoothed_values(paths, k) = to_host(step_payoffs(laws, Val(true), S_value, paths, k, args, mcBaseData, vb_mc))
    df = exp(-v_value(r) * dt)
    scale = v_value(mcProcess.underlying.S0)
    # Column k of S is t_{k-1}. Each path is exercised at t_{exercise_idx}, with cash flow (at that date) cash_flow.
    exercise_idx = fill(Nstep, Nsim)
    cash_flow = smoothed_values(1:Nsim, Nstep)
    for k = (Nstep-1):-1:1
        itm = findall(i -> v_value(FinancialMonteCarlo.payout(S_value_host[i, k], am_opt)) > 0, 1:Nsim)
        isempty(itm) && continue
        exercise_value = df .* smoothed_values(itm, k)
        continuation = @. cash_flow[itm] * df^(exercise_idx[itm] - k + 1)
        x = S_value_host[itm, k] ./ scale
        # All the paths start from S0: the continuation value is the mean.
        fitted = k == 1 ? fill(mean(continuation), length(itm)) : continuation_value.(Ref(continuation_coefficients(x, continuation)), x)
        exercised = exercise_value .> fitted
        @views @. cash_flow[itm[exercised]] = exercise_value[exercised] / df
        @views @. exercise_idx[itm[exercised]] = k
    end
    # Cash flows with derivatives, one smoothed payoff per path at its exercise date.
    cash = step_payoffs(laws, Val(false), S, 1:Nsim, exercise_idx, args, mcBaseData, vb_mc)
    return mean(@. cash * exp(-r * dt * $(on_device(mcBaseData, mcProcess, exercise_idx))))
end

# The entries (paths[j], cols[j]) of the matrix M (cols a column or one per path), gathered by linear index on the device
# of M: no scalar indexing of a GPU array.
function gather_at(M, paths, cols, mcBaseData)
    index = @. paths + $(size(M, 1)) * (cols - 1)
    return M[on_device(mcBaseData, M, collect(index))]
end

# The laws of one time step of the log price, to smooth the payoff over the step after the exercise decision: step_payoffs
# gives E[payoff(S_{t+dt}) | S_t] of the paths `paths` at the columns `cols` of S, undiscounted, with the derivatives in
# the density weight (values only with Val(true)). Lévy processes: the FFT density of the step, the same for every path.
struct FrozenFFT{L, V}
    law::L
    law_value::V
end
step_law(laws::FrozenFFT, ::Val{false}) = laws.law
step_law(laws::FrozenFFT, ::Val{true}) = laws.law_value
step_payoffs(laws::FrozenFFT, values, S, paths, cols, args, mcBaseData, vb_mc) = aug_paths(step_law(laws, values), args.option, vb_mc, gather_at(S, paths, cols, mcBaseData), mcBaseData)
# Heston: Gaussian, with the variance of the path at the start of the step.
struct FrozenHeston{R, M, N, D}
    rd::R
    V::M
    V_value::N
    dt::D
end
step_variances(laws::FrozenHeston, ::Val{false}) = (laws.rd, laws.V)
step_variances(laws::FrozenHeston, ::Val{true}) = (v_value(laws.rd), laws.V_value)
function step_payoffs(laws::FrozenHeston, values, S, paths, cols, args, mcBaseData, vb_mc)
    rd, V = step_variances(laws, values)
    S_prev = gather_at(S, paths, cols, mcBaseData)
    V_prev = gather_at(V, paths, cols, mcBaseData)
    dt = laws.dt
    mu = @. log(S_prev) + (rd - V_prev / 2) * dt
    sigma = @. sqrt(V_prev * dt)
    return lrm_paths(mu, sigma, args, mcBaseData, vb_mc)
end

# Paths and laws of one step for the frozen method, on the device of the configuration.
function frozen_laws(mcProcess, rfCurve, mcBaseData, T)
    f = FinancialFFT.density(mcProcess, rfCurve, FinancialMonteCarlo.with_mode(FFTDensityInverter(), mcBaseData.mode), T / mcBaseData.Nstep)
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mcBaseData, T)
    return S, FrozenFFT(fft_step_law(f), value_law(f))
end
function frozen_laws(mcProcess::FinancialMonteCarlo.HestonProcess, rfCurve, mcBaseData, T)
    dt = T / mcBaseData.Nstep
    rd = rfCurve.r - FinancialMonteCarlo.dividend(mcProcess)
    (V, S) = simulate_with_vol(mcProcess, rfCurve, mcBaseData, T)
    return S, FrozenHeston(rd, V, v_value.(V), dt)
end
frozen_laws(::FinancialMonteCarlo.ItoProcess, rfCurve, mcBaseData, T) = throw(ArgumentError("Vibrato on American options supports BlackScholesProcess, HestonProcess and Lévy processes."))

function pricer(mcProcess::BlackScholesProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, method::Vibrato, am_opt::FinancialMonteCarlo.AmericanPayoff)
    mcBaseData, vb_mc = method.mc, method.vb_mc
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    T = am_opt.T
    Nstep = mcBaseData.Nstep
    dt = T / Nstep
    σ = mcProcess.σ
    d = FinancialMonteCarlo.dividend(mcProcess)
    S0 = mcProcess.underlying.S0
    # Pass 1: the exercise regions, from Longstaff-Schwartz on independent paths with values only.
    value_process = BlackScholesProcess(v_value(σ), Underlying(v_value(S0), v_value(d)))
    # The regressions index single paths: they run on the host, on values (with CudaMode the paths come to the host).
    S_policy = to_host(FinancialMonteCarlo.simulate(value_process, ZeroRate(v_value(r)), mcBaseData, T))
    regions = exercise_regions(S_policy, am_opt, exp(-v_value(r) * dt), v_value(S0), v_value(σ) * sqrt(T))
    # Pass 2: paths on the half steps. Column 2k is s_k = (k - 1/2) dt, the exercise date t_k = k dt is the midpoint
    # of [s_k, s_{k+1}] and given its ends log S_{t_k} is Gaussian (Brownian bridge).
    # same sampling, mode and random state (the default constructor, with every field), twice the steps
    mc_half = FinancialMonteCarlo.MonteCarloConfiguration(mcBaseData.Nsim, 2 * Nstep, mcBaseData.sampling, mcBaseData.mode, mcBaseData.rng, mcBaseData.reseed)
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mc_half, T)
    args = LRMInterfaceArguments(init_lrm_vec(vb_mc, mcBaseData), am_opt)
    drift = r - d - drift_correction(mcProcess)
    bridge_std = σ * sqrt(dt) / 2
    # Vibrato over the last half step for the payoff at maturity, then backwards over each exercise date, all the paths
    # at once (on the device of the paths in CudaMode).
    W = exp(-r * T) .* lrm_paths((@views @muladd @. log(S[:, 2*Nstep]) + drift * dt / 2), σ * sqrt(dt / 2), args, mcBaseData, vb_mc)
    for k = (Nstep-1):-1:1
        m = @views @. (log(S[:, 2k]) + log(S[:, 2k+2])) / 2
        W = bridge_step(W, m, bridge_std, regions[k], am_opt, exp(-r * dt * k), mcBaseData)
    end
    return mean(W)
end

# Longstaff-Schwartz backward induction on values only. Returns, for each exercise date t_k, the exercise region
# {payoff > 0, payoff > continuation value} as a list of intervals in log price.
function exercise_regions(S, am_opt, df, scale, width)
    Nstep = size(S, 2) - 1
    payoff(s) = v_value(FinancialMonteCarlo.payout(s, am_opt))
    cash_flow = payoff.(S[:, end])
    exercise_idx = fill(Nstep, size(S, 1))
    regions = Vector{Vector{Tuple{eltype(S), eltype(S)}}}(undef, Nstep - 1)
    for k = (Nstep-1):-1:1
        S_k = S[:, k+1]
        itm = findall(s -> payoff(s) > 0, S_k)
        regions[k] = Tuple{eltype(S), eltype(S)}[]
        isempty(itm) && continue
        continuation = @. cash_flow[itm] * df^(exercise_idx[itm] - k)
        α = continuation_coefficients(S_k[itm] ./ scale, continuation)
        exercised = [payoff(s) > continuation_value(α, s / scale) for s in S_k[itm]]
        @views @. cash_flow[itm[exercised]] = payoff(S_k[itm[exercised]])
        @views @. exercise_idx[itm[exercised]] = k
        exercise_gain(x) = payoff(exp(x)) > 0 ? payoff(exp(x)) - continuation_value(α, exp(x) / scale) : -one(x)
        regions[k] = positive_intervals(exercise_gain, log(scale) - 12width, log(scale) + 12width)
    end
    return regions
end

# Intervals where f > 0 on a grid of [lo, hi], with the ends refined by bisection; intervals touching lo or hi are unbounded.
function positive_intervals(f, lo, hi; npts = 4001)
    xs = range(lo, hi, length = npts)
    inside = f.(xs) .> 0
    function bisect(a, b)
        for _ = 1:80
            m = (a + b) / 2
            (f(m) > 0) == (f(a) > 0) ? (a = m) : (b = m)
        end
        return (a + b) / 2
    end
    intervals = Tuple{eltype(xs), eltype(xs)}[]
    j = 1
    while j <= npts
        if inside[j]
            first_in = j
            while j < npts && inside[j+1]
                j += 1
            end
            push!(intervals, (first_in == 1 ? -Inf : bisect(xs[first_in-1], xs[first_in]), j == npts ? Inf : bisect(xs[j], xs[j+1])))
        end
        j += 1
    end
    return intervals
end

const BRIDGE_NODES, BRIDGE_WEIGHTS = gausslegendre(48)
const BRIDGE_Z = GAUSS_LEGENDRE_Z

# ∫_{za}^{zb} h(x) q(x; m, s) dx for the Gaussian density q, with x = v_value(m) + v_value(s) z: the points carry no
# derivative, all the derivatives are in the density weight, as in lrm_interface_vec!.
function lrm_segment(za, zb, m, s, h)
    c = (za + zb) / 2
    hw = (zb - za) / 2
    return sum(BRIDGE_WEIGHTS[j] * hw * bridge_integrand(c + hw * BRIDGE_NODES[j], m, s, h) for j in eachindex(BRIDGE_NODES))
end
function bridge_integrand(z, m, s, h)
    x = @muladd v_value(m) + v_value(s) * z
    return exp(log_density_default(z)) * h(x) * v_mod(log_density(x, m, s))
end

# Value at the previous date of a path that, at the exercise date (bridge mean m, standard deviation s), is exercised
# in region (discounted payoff) and otherwise continues with value W.
# The integrands of the bridge as callable structs, and the sum over the segments with sum(f, segments): a generator
# whose elements build a closure is not inferred on Julia 1.10.
struct PayoffAt{O}
    option::O
end
(p::PayoffAt)(x) = payout_f(x, p.option)
struct UnitPayoff end
(::UnitPayoff)(x) = 1
segments_integral(segments, m, s, h) = sum(segment -> lrm_segment(segment[1], segment[2], m, s, h), segments)

function bridge_vibrato(W, m, s, region, am_opt, discount)
    mv = v_value(m)
    sv = v_value(s)
    exercise = Tuple{typeof(mv), typeof(mv)}[]
    for (lo, hi) in region
        a = max((lo - mv) / sv, -BRIDGE_Z)
        b = min((hi - mv) / sv, BRIDGE_Z)
        a < b && push!(exercise, (a, b))
    end
    # All the bridge mass in the continuation region.
    isempty(exercise) && return W
    payoff_integral = segments_integral(exercise, m, s, PayoffAt(am_opt))
    continuation = Tuple{typeof(mv), typeof(mv)}[]
    a = -BRIDGE_Z
    for (lo, hi) in exercise
        lo > a && push!(continuation, (a, lo))
        a = hi
    end
    a < BRIDGE_Z && push!(continuation, (a, BRIDGE_Z))
    stay = isempty(continuation) ? zero(payoff_integral) : segments_integral(continuation, m, s, UnitPayoff())
    return discount * payoff_integral + W * stay
end

# One exercise date for all the paths: bridge_vibrato path by path on the host.
bridge_step(W, m, s, region, am_opt, discount, mcBaseData) = bridge_vibrato.(W, m, s, Ref(region), Ref(am_opt), discount)
# On the GPU: the exercise region is the same for every path, so its pieces are the same in number on every path (each
# interval of the region, clamped to ±BRIDGE_Z in standard deviations of the path, and the gaps between them); an empty
# piece has zero width and adds zero, as in device_option_paths. Each piece is one broadcast of paths × nodes.
function bridge_step(W, m, s, region, am_opt, discount, mcBaseData::CudaConfig)
    isempty(region) && return W
    x_row = node_row(m, BRIDGE_NODES)
    w_row = node_row(m, BRIDGE_WEIGHTS)
    mv = v_value.(m)
    sv = v_value(s)
    cut(y) = @. clamp((y - mv) / sv, -BRIDGE_Z, BRIDGE_Z)
    starts = [cut(lo) for (lo, _) in region]
    ends = [cut(hi) for (_, hi) in region]
    segment(za, zb, h) = device_bridge_segment(za, zb, m, s, h, x_row, w_row)
    payoff_integral = sum(segment(za, zb, PayoffAt(am_opt)) for (za, zb) in zip(starts, ends))
    # the continuation region: the gaps [-Z, a_1], [b_1, a_2], ..., [b_n, Z]
    lower = fill!(similar(mv), -BRIDGE_Z)
    upper = fill!(similar(mv), BRIDGE_Z)
    stay = sum(segment(za, zb, UnitPayoff()) for (za, zb) in zip([[lower]; ends], [starts; [upper]]))
    # a path whose bridge mass is all in the continuation region keeps W, as in bridge_vibrato
    exercised = reduce((e, (za, zb)) -> e .| (za .< zb), zip(starts, ends); init = false)
    return @. bridge_value(exercised, W, payoff_integral, stay, discount)
end
device_bridge_segment(za, zb, m, s, h, x_row, w_row) = (h_ref = Ref(h); sum_nodes(@muladd @. w_row * (zb - za) / 2 * bridge_integrand((za + zb) / 2 + (zb - za) / 2 * x_row, m, s, h_ref)))
function bridge_value(exercised, W, payoff_integral, stay, discount)
    value = @muladd discount * payoff_integral + W * stay
    return ifelse(exercised, value, oftype(value, W))
end
