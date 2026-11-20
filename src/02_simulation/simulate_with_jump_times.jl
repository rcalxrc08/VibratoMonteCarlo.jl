using FinancialMonteCarlo, Random, Distributions, Sobol, StatsFuns;

# The last jump of each path of Kou and Merton, for VibratoSaltando and ConditionalSaltando. The diffusion is simulated by
# FinancialMonteCarlo where the paths live (CURAND on the GPU), the jumps are drawn by FinancialMonteCarlo into a JumpList
# (path, first grid point, size) with its schemes for every sampling (add_jumps!: exponential spacings for plain MC, the
# uniforms of JumpUniforms, jump_count and jump_quantile for AntitheticMC and SobolMC), and the paths never leave the
# device: the last jump of a path only needs the diffusion at one grid point and the sum of its jumps, gathered there.
#
# A jump in step k is applied at grid point k + 1, as in FinancialMonteCarlo. The conditioning is on the grid point
# before the step containing the last jump, plus every jump of the path (those of that step included): the diffusion
# left to T is then exactly Gaussian and lasts at least one step. (The grid point of the jump already contains diffusion
# beyond the jump time: using it with that time counts that diffusion twice.)

# Plain MC and AntitheticMC: FinancialMonteCarlo's own schemes, with a JumpList as the sink. Antithetic pairs reflect the
# uniforms of the jumps (1 - u) on the path i + Nsim / 2, as FinancialMonteCarlo pairs the normals.
draw_jumps!(jumps, mcProcess, mcBaseData, T) = FinancialMonteCarlo.add_jumps!(jumps, mcProcess, mcBaseData, T)
# Quasi Monte Carlo: the jumps take the coordinates Nstep + 1, ..., Nstep + 1 + 2M of the Sobol point of the path, whose
# first Nstep coordinates are the ones FinancialMonteCarlo used for its diffusion (same dimensions, same skip): each path
# is one point of a Sobol sequence of dimension Nstep + 1 + 2M, as in FinancialMonteCarlo's Sobol Kou and Merton.
function draw_jumps!(jumps, mcProcess, mcBaseData::FinancialMonteCarlo.SobolConfig, T)
    J = FinancialMonteCarlo.JumpUniforms(v_value(mcProcess.λ * T), Float64)
    Nstep = mcBaseData.Nstep
    dimension = Nstep + FinancialMonteCarlo.jump_uniforms(J)
    seq = FinancialMonteCarlo.sobol_skip!(SobolSeq(dimension), Nstep * mcBaseData.Nsim)
    point = Vector{Float64}(undef, dimension)
    u = view(point, (Nstep+1):dimension)
    for i = 1:mcBaseData.Nsim
        next!(seq, point)
        FinancialMonteCarlo.add_jumps_from_uniforms!(jumps, i, u, J, mcProcess, mcBaseData)
    end
    return jumps
end

# The grid point where the last jump of each path is applied (2 without jumps: the start of the path) and the sum of its
# jumps, in one pass over the list (not over the paths).
function last_jumps(jumps::FinancialMonteCarlo.JumpList, Nsim, zero_drift)
    last_step = fill(Int32(2), Nsim)
    jump_sum = fill(zero_drift, Nsim)
    @inbounds for k in eachindex(jumps.path)
        i = jumps.path[k]
        last_step[i] = max(last_step[i], jumps.step[k])
        jump_sum[i] += jumps.size[k]
    end
    return last_step, jump_sum
end

function simulate_log_with_jump_times(mcProcess::finiteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess <: FinancialMonteCarlo.FiniteActivityProcess, numb <: Number}
    r = rfCurve.r
    d = FinancialMonteCarlo.dividend(mcProcess)
    σ = mcProcess.σ
    λ = mcProcess.λ
    # r - d - ψ(-i)
    drift_rn = (r - d) - drift_correction(mcProcess)
    zero_drift = drift_rn * 0 * T
    Nsim = mcBaseData.Nsim
    dt = T / mcBaseData.Nstep
    X = FinancialMonteCarlo.get_matrix_type(mcBaseData, mcProcess, typeof(zero_drift))
    FinancialMonteCarlo.simulate!(X, BrownianMotion(σ, drift_rn), mcBaseData, T)
    # room for the expected number of jumps and 4 standard deviations, as FinancialMonteCarlo
    expected = Nsim * v_value(λ * T)
    jumps = FinancialMonteCarlo.JumpList{typeof(zero_drift)}(ceil(Int, expected + 4 * sqrt(expected)))
    draw_jumps!(jumps, mcProcess, mcBaseData, T)
    last_step, jump_sum = last_jumps(jumps, Nsim, zero_drift)
    # the diffusion at the grid point before the step of the last jump, X[i, last_step[i] - 1], gathered by linear index
    grid_point = on_device(mcBaseData, mcProcess, @. (1:Nsim) + Nsim * (last_step - 2))
    zero_time = 0 * λ
    T_jumps = on_device(mcBaseData, mcProcess, @. zero_time + (last_step - 2) * dt)
    X_f = @. log(mcProcess.underlying.S0) + $(X[grid_point]) + $(on_device(mcBaseData, mcProcess, jump_sum))
    return T_jumps, X_f
end
function simulate_with_jump_times(mcProcess::finiteActivityProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess <: FinancialMonteCarlo.FiniteActivityProcess, numb <: Number}
    T_jumps, X_f = simulate_log_with_jump_times(mcProcess, rfCurve, mcBaseData, T)
    S_jumps = @. exp(X_f)
    return T_jumps, S_jumps
end
