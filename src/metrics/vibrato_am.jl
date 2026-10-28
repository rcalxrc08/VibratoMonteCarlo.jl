using LinearAlgebra

# Continuation value regressed on 1, x, x^2 (Longstaff-Schwartz), values only.
continuation_coefficients(x, y) = hcat(one.(x), x, x .^ 2) \ y
continuation_value(α, x) = α[1] + x * (α[2] + x * α[3])

"""
    vibrato_am(mcProcess, rfCurve, mcBaseData, am_opt, vb_mc)

Price of an American option exercisable at the `Nstep` dates of the time grid, with its sensitivities.

The exercise policy comes from Longstaff-Schwartz (regression on `1, S, S^2` of the in the money paths, values only), so
it carries no derivative. How the derivatives get through the exercise decisions depends on the process:
- `BlackScholesProcess`: vibrato over every exercise date. The policy is fixed in price space and each path is simulated on
  the half steps, so that every exercise date is the midpoint of a Brownian bridge: the exercise indicator and the payoff
  are integrated against the bridge density with the likelihood ratio method. Unbiased for the policy at every order,
  including the second order sensitivities w.r.t. `S0`; the variance of the second order ones grows with `Nstep`.
- Lévy processes (FFT density of one step, requires `VibratoMonteCarloAnalytic`) and `HestonProcess` (Gaussian step with
  the simulated variance, `ZeroRate` and serial or antithetic configurations only): the decision to exercise at `t_k` is
  taken at `t_{k-1}`, the payoff is replaced by its vibrato expectation over the step and the decisions are frozen when
  differentiating (Pagès, Pironneau, Sall). Value and first order sensitivities are accurate, second order ones miss
  the effect of the moving exercise boundary and are biased.
"""
function vibrato_am(mcProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, am_opt::FinancialMonteCarlo.AmericanPayoff, vb_mc::AbstractVibrato)
    FinancialMonteCarlo.set_seed!(mcBaseData)
    r = rfCurve.r
    T = am_opt.T
    Nstep = mcBaseData.Nstep
    dt = T / Nstep
    S, law_at, value_law_at = frozen_laws(mcProcess, rfCurve, mcBaseData, T)
    S_value = v_value.(S)
    Nsim = size(S, 1)
    args = LRMInterfaceArguments(init_lrm_vec(vb_mc, mcBaseData), am_opt)
    smoothed_value(i, k) = smoothed_payoff(value_law_at(i, k), S_value[i, k], args, mcProcess, mcBaseData, vb_mc)
    df = exp(-v_value(r) * dt)
    scale = v_value(mcProcess.underlying.S0)
    # Column k of S is t_{k-1}. Each path is exercised at t_{exercise_idx}, with cash flow (at that date) cash_flow.
    exercise_idx = fill(Nstep, Nsim)
    cash_flow = [smoothed_value(i, Nstep) for i = 1:Nsim]
    for k = (Nstep-1):-1:1
        itm = findall(i -> v_value(FinancialMonteCarlo.payout(S_value[i, k], am_opt)) > 0, 1:Nsim)
        isempty(itm) && continue
        exercise_value = [df * smoothed_value(i, k) for i in itm]
        continuation = @. cash_flow[itm] * df^(exercise_idx[itm] - k + 1)
        x = S_value[itm, k] ./ scale
        # All the paths start from S0: the continuation value is the mean.
        fitted = k == 1 ? fill(mean(continuation), length(itm)) : continuation_value.(Ref(continuation_coefficients(x, continuation)), x)
        exercised = exercise_value .> fitted
        @views @. cash_flow[itm[exercised]] = exercise_value[exercised] / df
        @views @. exercise_idx[itm[exercised]] = k
    end
    # Cash flows with derivatives, one smoothed payoff per path.
    return mean(smoothed_payoff(law_at(i, exercise_idx[i]), S[i, exercise_idx[i]], args, mcProcess, mcBaseData, vb_mc) * exp(-r * dt * exercise_idx[i]) for i = 1:Nsim)
end

# Law of one time step of the log price, used to smooth the payoff over the step after the exercise decision:
# smoothed_payoff(law, S_prev, ...) = E[payoff(S_{t+dt}) | S_t = S_prev], undiscounted, with the derivatives in the density weight.
struct GaussianStep{D, S}
    drift::D
    sigma::S
end

smoothed_payoff(law::GaussianStep, S_prev, args, mcProcess, mcBaseData, vb_mc) = lrm_interface_vec!(log(S_prev) + law.drift, law.sigma, args, mcBaseData, vb_mc)
smoothed_payoff(law::DensitySpline, S_prev, args, mcProcess, mcBaseData, vb_mc) = lrm_interface_aug!(mcProcess, law, args.option, mcBaseData, vb_mc, S_prev)

# Paths and laws of one step for the frozen method: returns S, law_at(i, k) and its values only version, for path i at column k.
function frozen_laws(mcProcess, rfCurve, mcBaseData, T)
    (x, density_val) = density_unadapted(mcProcess, T / mcBaseData.Nstep, rfCurve.r, 18, 18.0)
    law = density_spline(x, density_val)
    law_value = density_spline(x, v_value.(density_val))
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mcBaseData, T)
    return S, (_, _) -> law, (_, _) -> law_value
end
function frozen_laws(mcProcess::FinancialMonteCarlo.HestonProcess, rfCurve, mcBaseData, T)
    dt = T / mcBaseData.Nstep
    rd = rfCurve.r - FinancialMonteCarlo.dividend(mcProcess)
    (V, S) = simulate_with_vol(mcProcess, rfCurve, mcBaseData, T)
    V_value = v_value.(V)
    law_at(i, k) = GaussianStep((rd - V[i, k] / 2) * dt, sqrt(V[i, k] * dt))
    value_law_at(i, k) = GaussianStep((v_value(rd) - V_value[i, k] / 2) * dt, sqrt(V_value[i, k] * dt))
    return S, law_at, value_law_at
end
frozen_laws(::FinancialMonteCarlo.ItoProcess, rfCurve, mcBaseData, T) = throw(ArgumentError("vibrato_am supports BlackScholesProcess, HestonProcess and Lévy processes."))

function vibrato_am(mcProcess::BlackScholesProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, am_opt::FinancialMonteCarlo.AmericanPayoff, vb_mc::AbstractVibrato)
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
    S_policy = FinancialMonteCarlo.simulate(value_process, ZeroRate(v_value(r)), mcBaseData, T)
    regions = exercise_regions(S_policy, am_opt, exp(-v_value(r) * dt), v_value(S0), v_value(σ) * sqrt(T))
    # Pass 2: paths on the half steps. Column 2k is s_k = (k - 1/2) dt, the exercise date t_k = k dt is the midpoint
    # of [s_k, s_{k+1}] and given its ends log S_{t_k} is Gaussian (Brownian bridge).
    mc_half = MonteCarloConfiguration(mcBaseData.Nsim, 2 * Nstep, mcBaseData.monteCarloMethod, mcBaseData.parallelMode, mcBaseData.init_rng)
    S = FinancialMonteCarlo.simulate(mcProcess, rfCurve, mc_half, T)
    args = LRMInterfaceArguments(init_lrm_vec(vb_mc, mcBaseData), am_opt)
    drift = r - d - σ^2 / 2
    bridge_std = σ * sqrt(dt) / 2
    function path_value(i)
        # Vibrato over the last half step for the payoff at maturity, then backwards over each exercise date.
        W = exp(-r * T) * lrm_interface_vec!(log(S[i, 2*Nstep]) + drift * dt / 2, σ * sqrt(dt / 2), args, mcBaseData, vb_mc)
        for k = (Nstep-1):-1:1
            W = bridge_vibrato(W, (log(S[i, 2k]) + log(S[i, 2k+2])) / 2, bridge_std, regions[k], am_opt, exp(-r * dt * k))
        end
        return W
    end
    return mean(path_value(i) for i = 1:size(S, 1))
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

# Gauss-Legendre nodes and weights on [-1, 1].
function gauss_legendre(n)
    x = zeros(n)
    w = zeros(n)
    for i = 1:n
        z = cos(π * (i - 0.25) / (n + 0.5))
        dp = 1.0
        for _ = 1:100
            p1, p2 = 1.0, 0.0
            for j = 1:n
                p1, p2 = ((2j - 1) * z * p1 - (j - 1) * p2) / j, p1
            end
            dp = n * (z * p1 - p2) / (z^2 - 1)
            Δ = p1 / dp
            z -= Δ
            abs(Δ) < 1e-15 && break
        end
        x[i] = z
        w[i] = 2 / ((1 - z^2) * dp^2)
    end
    return x, w
end
const BRIDGE_NODES, BRIDGE_WEIGHTS = gauss_legendre(48)
# The bridge mass beyond ±BRIDGE_Z standard deviations is below 1e-16.
const BRIDGE_Z = 8.5

# ∫_{za}^{zb} h(x) q(x; m, s) dx for the Gaussian density q, with x = v_value(m) + v_value(s) z: the points carry no
# derivative, all the derivatives are in the density weight, as in lrm_interface_vec!.
function lrm_segment(za, zb, m, s, h)
    c = (za + zb) / 2
    hw = (zb - za) / 2
    return sum(begin
        z = c + hw * BRIDGE_NODES[j]
        x = v_value(m) + v_value(s) * z
        BRIDGE_WEIGHTS[j] * hw * exp(log_density_default(z)) * h(x) * v_mod(log_density(x, m, s))
    end for j in eachindex(BRIDGE_NODES))
end

# Value at the previous date of a path that, at the exercise date (bridge mean m, standard deviation s), is exercised
# in region (discounted payoff) and otherwise continues with value W.
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
    payoff_integral = sum(lrm_segment(a, b, m, s, x -> payout_f(x, am_opt)) for (a, b) in exercise)
    continuation = Tuple{typeof(mv), typeof(mv)}[]
    a = -BRIDGE_Z
    for (lo, hi) in exercise
        lo > a && push!(continuation, (a, lo))
        a = hi
    end
    a < BRIDGE_Z && push!(continuation, (a, BRIDGE_Z))
    stay = isempty(continuation) ? zero(payoff_integral) : sum(lrm_segment(a, b, m, s, _ -> 1) for (a, b) in continuation)
    return discount * payoff_integral + W * stay
end
