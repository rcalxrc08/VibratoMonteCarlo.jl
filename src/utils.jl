abstract type AbstractVibrato end
using FinancialToolbox
"""
    VibratoMonteCarloStandard(Nsim_lrm = 200)

Inner integration by Monte Carlo: `Nsim_lrm` standard normal samples for every outer path, drawn again for each path from
the generator of the `MonteCarloConfiguration`. With `AntitheticMC` they come in antithetic pairs, with `SobolMC`
they are a one dimensional Sobol sequence randomly shifted for every path (Cranley-Patterson). Gaussian inner laws only.
"""
struct VibratoMonteCarloStandard <: AbstractVibrato
    Nsim_lrm::Int
    function VibratoMonteCarloStandard(Nsim_lrm::Int = 200)
        @assert Nsim_lrm > 0 "Error, number of simulation for VibratoMonteCarloStandard must be positive."
        return new(Nsim_lrm)
    end
end

function init_lrm_vec(x::VibratoMonteCarloStandard, mcBaseData)
    Z = Array{Float64}(undef, x.Nsim_lrm)
    randn!(FinancialMonteCarlo.host_rng(mcBaseData), Z)
    return Z
end

"""
    VibratoMonteCarloAnalytic(Npoints = 200, x_min = -7.0, x_max = 7.0)

Inner integration by quadrature on a uniform grid of `Npoints` in `[x_min, x_max]`, in standard deviations of the log
increment around its mean. The only inner integration of the FFT density (Lévy processes with `Vibrato` and `PureLRM`).
"""
struct VibratoMonteCarloAnalytic <: AbstractVibrato
    Npoints::Int
    x_min::Float64
    x_max::Float64
    # ±5 standard deviations miss part of the score terms in the upper tail, where a call payoff is large: volga was 3-45% off on VG.
    function VibratoMonteCarloAnalytic(Npoints::Int = 200, x_min::Float64 = -7.0, x_max::Float64 = 7.0)
        @assert Npoints > 0 "Error, number of integration points for VibratoMonteCarloAnalytic must be positive."
        return new(Npoints, x_min, x_max)
    end
end
Base.broadcastable(x::T) where {T <: AbstractVibrato} = Ref(x)

function init_lrm_vec(x::VibratoMonteCarloAnalytic, ::Any)
    return range(x.x_min, length = x.Npoints, stop = x.x_max)
end

"""
    VibratoMonteCarlo.VibratoMonteCarloGaussHermite(Npoints = 200)

Inner integration by Gauss-Hermite quadrature on `Npoints` nodes, Gaussian inner laws only. It converges slowly across the
kink or the jump of a payoff: prefer `VibratoMonteCarloGaussLegendre`.
"""
struct VibratoMonteCarloGaussHermite <: AbstractVibrato
    Npoints::Int
    function VibratoMonteCarloGaussHermite(Npoints::Int = 200)
        @assert Npoints > 0 "Error, number of integration points for VibratoMonteCarloGaussHermite must be positive."
        return new(Npoints)
    end
end

init_lrm_vec(x::VibratoMonteCarloGaussHermite, ::Any) = gausshermite(x.Npoints)

# Gauss-Legendre in the standardized variable, on [-GAUSS_LEGENDRE_Z, GAUSS_LEGENDRE_Z + sigma] split at the payoff singular points
# (singular_points), with Npoints nodes per piece. Every piece has a smooth integrand, so the rule converges exponentially:
# Gauss-Hermite instead converges slowly across the kink or the jump of the payoff.
"""
    VibratoMonteCarloGaussLegendre(Npoints = 48)

Inner integration by Gauss-Legendre quadrature on `[-8.5, 8.5 + σ]` standard deviations, split at the singular points of
the payoff (the strikes) with `Npoints` nodes on each piece. Every piece has a smooth integrand, so the rule converges
exponentially: about 1e-14 per path. Gaussian inner laws only (Black and Scholes, Heston, `VibratoSaltando`, `PureLRM`
on Black and Scholes, American options on Black and Scholes).
"""
struct VibratoMonteCarloGaussLegendre <: AbstractVibrato
    Npoints::Int
    function VibratoMonteCarloGaussLegendre(Npoints::Int = 48)
        @assert Npoints > 0 "Error, number of integration points for VibratoMonteCarloGaussLegendre must be positive."
        return new(Npoints)
    end
end
# The Gaussian mass beyond ±GAUSS_LEGENDRE_Z standard deviations is below 1e-16.
const GAUSS_LEGENDRE_Z = 8.5

# Nodes and weights on [-1, 1].
init_lrm_vec(x::VibratoMonteCarloGaussLegendre, ::Any) = gausslegendre(x.Npoints)

get_unique_maturity(euopt) = euopt.T

# GPU (FinancialMonteCarlo.CudaMode, with `using CUDA`): FinancialMonteCarlo simulates on the device, the inner integration
# runs there too (lrm_paths, aug_paths), and only the regressions of American options take the values of the paths to
# the host. Every transfer is one copy of a whole array, so no GPU array is indexed by scalars. AbstractCudaMode exists
# without CUDA, so no extension is needed.
const CudaConfig = FinancialMonteCarlo.MonteCarloConfiguration{<:FinancialMonteCarlo.AbstractSampling, <:FinancialMonteCarlo.AbstractCudaMode}
to_host(x::Array) = x
to_host(x::AbstractArray) = Array(x)
# x on the device of the mode of the configuration (CuArray with the element type kept, or Array).
on_device(mcBaseData, mcProcess, x) = FinancialMonteCarlo.to_mode(x, mcBaseData.mode)
# Standard normals in place: from the host generator, or from CURAND for a CuArray (FinancialMonteCarlo's CUDA extension,
# seeded by set_seed!).
normals!(mcBaseData, e) = FinancialMonteCarlo.draw_normals!(e, mcBaseData)
# Mean over the paths of lrm_paths: a vector (of vectors for a smile on the CPU), or a matrix with one column per option.
path_mean(x::AbstractVector) = mean(x)
path_mean(x::AbstractMatrix) = vec(to_host(mean(x; dims = 1)))

# log E[exp(X_1)] of the driftless log price, from the characteristic exponent of FinancialFFT: the risk neutral drift of
# log S is r - d - drift_correction (σ² / 2 for Black-Scholes, ψ for VG and NIG, σ² / 2 + λ E[e^J - 1] with jumps).
drift_correction(mcProcess) = FinancialFFT.characteristic_exponent_i(1, mcProcess)

function analytic_pricer(S_jump, eu_opt::EuropeanOption, r, T, T_jump, σ, d, drift_rn)
    blsprice(S_jump, eu_opt.K, r, T - T_jump, σ, d + drift_rn)
end
function analytic_pricer(S_jump, eu_opt::BinaryEuropeanOption, r, T, T_jump, σ, d, drift_rn)
    blsbin(S_jump, eu_opt.K, r, T - T_jump, σ, d + drift_rn)
end
