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
    function VibratoMonteCarloAnalytic(Npoints::Int = 200, x_min::Float64 = -5.0, x_max::Float64 = 5.0)
        @assert Npoints > 0 "Error, number of integration points for VibratoMonteCarloAnalytic must be positive."
        return new(Npoints, x_min, x_max)
    end
end
Base.broadcastable(x::T) where {T <: AbstractVibrato} = Ref(x)

function init_lrm_vec(x::VibratoMonteCarloAnalytic, ::Any)
    return range(x.x_min, length = x.Npoints, stop = x.x_max)
end

get_unique_maturity(euopt) = euopt.T
get_unique_maturity(euopts::Array) = only(unique([get_unique_maturity(euopt) for euopt in euopts]))

function analytic_pricer(S_jump, eu_opt::EuropeanOption, r, T, T_jump, σ, d, drift_rn)
    blsprice(S_jump, eu_opt.K, r, T - T_jump, σ, d + drift_rn)
end
function analytic_pricer(S_jump, eu_opt::BinaryEuropeanOption, r, T, T_jump, σ, d, drift_rn)
    blsbin(S_jump, eu_opt.K, r, T - T_jump, σ, d + drift_rn)
end