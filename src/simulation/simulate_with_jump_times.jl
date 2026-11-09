using FinancialMonteCarlo, Random, Distributions, Sobol, StatsFuns;

function generate_last_jump_time_and_value(row_X_i, λ, zero_time, zero_drift, T, mcProcess, dt, mcBaseData)
    t_i = randexp(FinancialMonteCarlo.host_rng(mcBaseData)) / λ
    t_i_1 = t_i / T
    T_jump = zero_time
    x_jump = zero_drift
    last_idx = 0
    jumps_in_step = zero_drift
    while t_i_1 < 1
        jump_size = FinancialMonteCarlo.compute_jump_size(mcProcess, mcBaseData)
        jump_idx = ceil(UInt32, t_i / dt) + true
        @views @. row_X_i[jump_idx:end] += jump_size #add jump component
        jumps_in_step = jump_idx == last_idx ? jumps_in_step + jump_size : jump_size
        last_idx = jump_idx
        # Condition on the grid point before the step containing the last jump, plus the jumps of that step:
        # the diffusion left to T is then exactly Gaussian and lasts at least one step.
        # (The value at jump_idx already contains diffusion beyond t_i, using it with t_i counts that diffusion twice.)
        T_jump = (jump_idx - 2) * dt
        x_jump = row_X_i[jump_idx-1] + jumps_in_step
        t_i += randexp(FinancialMonteCarlo.host_rng(mcBaseData)) / λ
        t_i_1 = t_i / T
    end
    return T_jump, x_jump
end

# Plain Monte Carlo: the jump times are exponential spacings and the sizes come from FinancialMonteCarlo, both from the rng.
function add_last_jumps!(T_jumps, S_jumps, X, λ, zero_time, zero_drift, T, mcProcess, dt, mcBaseData)
    @inbounds for i = 1:mcBaseData.Nsim
        T_jumps[i], S_jumps[i] = generate_last_jump_time_and_value(view(X, i, :), λ, zero_time, zero_drift, T, mcProcess, dt, mcBaseData)
    end
    return
end

# AntitheticMC and SobolMC need the jumps of a path from a fixed number of uniforms, while the number of jumps is random
# (#30): the count N_T is the Poisson quantile at u[1], the jump times are the uniforms u[2:N_T+1], sorted (given N_T they
# are the order statistics of N_T uniforms on [0, T]), and the sizes the quantiles of the jump law at u[M+2:M+N_T+1].
# M is the Poisson quantile at 1 - 1e-9: a path with more jumps draws all of them from the rng, as plain MC.
# The Poisson distribution function up to M is tabulated once, and the times are sorted in a buffer: no allocation per
# path (the quantile of Distributions and the vectors per path made the simulation 3-5 times slower).
struct JumpUniforms{L}
    law::L
    M::Int
    cdf::Vector{Float64}
    times::Vector{Float64}
end
function JumpUniforms(λT)
    law = Poisson(λT)
    M = max(1, quantile(law, 1 - 1e-9))
    return JumpUniforms(law, M, cdf.(law, 0:M), Vector{Float64}(undef, M))
end
jump_uniforms(J::JumpUniforms) = 1 + 2 * J.M
# The smallest N with P(N_T <= N) >= u, as quantile(law, u).
jump_count(J::JumpUniforms, u) = u <= last(J.cdf) ? searchsortedfirst(J.cdf, u) - 1 : quantile(J.law, u)
# Inversion of the jump laws of FinancialMonteCarlo's compute_jump_size: normal for Merton, double exponential for Kou
# (downwards with probability 1 - p and rate λ₋, upwards with rate λ₊). Kou's quantile is continuous at u = 1 - p.
jump_quantile(mcProcess::MertonProcess, u) = @muladd mcProcess.μ_jump + mcProcess.σ_jump * norminvcdf(u)
function jump_quantile(mcProcess::KouProcess, u)
    p = mcProcess.p
    q = 1 - p
    return u < v_value(q) ? log(u / q) / mcProcess.λ₋ : -log((1 - u) / p) / mcProcess.λ₊
end
# A uniform of a Sobol point or of a reflection can be 0 or 1, where the quantiles are infinite.
open_unit(u) = clamp(u, eps(), 1 - eps())

# The jumps of one path from its uniforms u, added to the row of its diffusion: the same last jump as
# generate_last_jump_time_and_value. The time of a jump only fixes its grid step: a jump in step k is applied at grid
# point k + 1, as in FinancialMonteCarlo. The sizes are independent of the times, so the k-th smallest time takes the
# k-th size uniform.
function last_jump_from_uniforms(row_X_i, u, J::JumpUniforms, zero_time, zero_drift, mcProcess, dt, mcBaseData)
    N = jump_count(J, u[1])
    N > J.M && return last_jump_from_rng(row_X_i, N, zero_time, zero_drift, mcProcess, dt, mcBaseData)
    times = view(J.times, 1:N)
    copyto!(times, view(u, 2:(N+1)))
    sort!(times)
    return add_sorted_jumps!(row_X_i, times, k -> jump_quantile(mcProcess, open_unit(u[J.M+1+k])), zero_time, zero_drift, dt)
end
# More than M jumps (probability 1e-9): times and sizes from the rng.
function last_jump_from_rng(row_X_i, N, zero_time, zero_drift, mcProcess, dt, mcBaseData)
    times = sort!(rand(FinancialMonteCarlo.host_rng(mcBaseData), N))
    return add_sorted_jumps!(row_X_i, times, _ -> FinancialMonteCarlo.compute_jump_size(mcProcess, mcBaseData), zero_time, zero_drift, dt)
end
function add_sorted_jumps!(row_X_i, times, size_of, zero_time, zero_drift, dt)
    Nstep = length(row_X_i) - 1
    T_jump = zero_time
    x_jump = zero_drift
    last_idx = 0
    jumps_in_step = zero_drift
    for (k, t) in enumerate(times)
        jump_idx = max(ceil(Int, t * Nstep), 1) + 1
        jump_size = size_of(k)
        @views @. row_X_i[jump_idx:end] += jump_size
        jumps_in_step = jump_idx == last_idx ? jumps_in_step + jump_size : jump_size
        last_idx = jump_idx
        T_jump = (jump_idx - 2) * dt
        x_jump = row_X_i[jump_idx-1] + jumps_in_step
    end
    return T_jump, x_jump
end

# Antithetic pairs as FinancialMonteCarlo's diffusion: the path i + Nsim / 2 has the normals of the path i with the
# opposite sign, and here the uniforms 1 - u of its jumps.
function add_last_jumps!(T_jumps, S_jumps, X, λ, zero_time, zero_drift, T, mcProcess, dt, mcBaseData::FinancialMonteCarlo.MonteCarloConfiguration{<:FinancialMonteCarlo.AntitheticMC})
    J = JumpUniforms(v_value(λ * T))
    u = Vector{Float64}(undef, jump_uniforms(J))
    u_anti = similar(u)
    half = div(mcBaseData.Nsim, 2)
    @inbounds for i = 1:half
        rand!(FinancialMonteCarlo.host_rng(mcBaseData), u)
        @. u_anti = 1 - u
        T_jumps[i], S_jumps[i] = last_jump_from_uniforms(view(X, i, :), u, J, zero_time, zero_drift, mcProcess, dt, mcBaseData)
        T_jumps[i+half], S_jumps[i+half] = last_jump_from_uniforms(view(X, i + half, :), u_anti, J, zero_time, zero_drift, mcProcess, dt, mcBaseData)
    end
    return
end

# Quasi Monte Carlo: the jumps take the coordinates Nstep + 1, ..., Nstep + 1 + 2M of the Sobol point of the path, whose
# first Nstep coordinates are the ones FinancialMonteCarlo used for its diffusion (same dimensions, same skip): each path
# is one point of a Sobol sequence of dimension Nstep + 1 + 2M.
function add_last_jumps!(T_jumps, S_jumps, X, λ, zero_time, zero_drift, T, mcProcess, dt, mcBaseData::FinancialMonteCarlo.SobolConfig)
    J = JumpUniforms(v_value(λ * T))
    Nstep = mcBaseData.Nstep
    dimension = Nstep + jump_uniforms(J)
    seq = FinancialMonteCarlo.sobol_skip!(SobolSeq(dimension), Nstep * mcBaseData.Nsim)
    point = Vector{Float64}(undef, dimension)
    u = view(point, (Nstep+1):dimension)
    @inbounds for i = 1:mcBaseData.Nsim
        next!(seq, point)
        T_jumps[i], S_jumps[i] = last_jump_from_uniforms(view(X, i, :), u, J, zero_time, zero_drift, mcProcess, dt, mcBaseData)
    end
    return
end

function simulate_log_with_jump_times(mcProcess::finiteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess <: FinancialMonteCarlo.FiniteActivityProcess, numb <: Number}
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    λ = mcProcess.λ

    ####Simulation
    ## Simulate
    # r-d-psi(-i)
    drift_rn = (r - d) - drift_correction(mcProcess)
    zero_drift = drift_rn * 0 * T
    Nsim = mcBaseData.Nsim
    Nsteps = mcBaseData.Nstep
    X = FinancialMonteCarlo.get_matrix_type(mcBaseData, mcProcess, typeof(zero_drift))
    FinancialMonteCarlo.simulate!(X, BrownianMotion(σ, drift_rn), mcBaseData, T)
    # The jumps are added path by path with the CPU rng: on the host, then the results go back to the device.
    X = to_host(X)
    S0 = mcProcess.underlying.S0
    dt = T / Nsteps
    T_jumps = Array{typeof(λ)}(undef, Nsim)
    S_jumps = Array{typeof(zero_drift)}(undef, Nsim)
    zero_time = 0 * λ
    add_last_jumps!(T_jumps, S_jumps, X, λ, zero_time, zero_drift, T, mcProcess, dt, mcBaseData)
    log_S0 = log(S0)
    X_f = @. log_S0 + S_jumps
    return on_device(mcBaseData, mcProcess, T_jumps), on_device(mcBaseData, mcProcess, X_f)
end
function simulate_with_jump_times(mcProcess::finiteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess <: FinancialMonteCarlo.FiniteActivityProcess, numb <: Number}
    T_jumps, X_f = simulate_log_with_jump_times(mcProcess, rfCurve, mcBaseData, T)
    S_jumps = @. exp(X_f)
    return T_jumps, S_jumps
end
