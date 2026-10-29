abstract type AbstractVibrato end
using FinancialToolbox
struct VibratoMonteCarloStandard <: AbstractVibrato
    Nsim_lrm::Int
    mode::FinancialMonteCarlo.SerialMode
    function VibratoMonteCarloStandard(Nsim_lrm::Int = 200, mode::FinancialMonteCarlo.SerialMode = FinancialMonteCarlo.SerialMode())
        @assert Nsim_lrm > 0 "Error, number of simulation for VibratoMonteCarloStandard must be positive."
        return new(Nsim_lrm, mode)
    end
end

function init_lrm_vec(x::VibratoMonteCarloStandard, mcBaseData)
    Z = Array{Float64}(undef, x.Nsim_lrm)
    randn!(mcBaseData.parallelMode.rng, Z)
    return Z
end

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
get_unique_maturity(euopts::Array) = only(unique([get_unique_maturity(euopt) for euopt in euopts]))

function analytic_pricer(S_jump, eu_opt::EuropeanOption, r, T, T_jump, σ, d, drift_rn)
    blsprice(S_jump, eu_opt.K, r, T - T_jump, σ, d + drift_rn)
end
function analytic_pricer(S_jump, eu_opt::BinaryEuropeanOption, r, T, T_jump, σ, d, drift_rn)
    blsbin(S_jump, eu_opt.K, r, T - T_jump, σ, d + drift_rn)
end
