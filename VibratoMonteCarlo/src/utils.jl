abstract type AbstractVibrato end

struct VibratoMonteCarloStandard <: AbstractVibrato
    Nsim_lrm::Int
    mode::FinancialMonteCarlo.SerialMode
    function VibratoMonteCarloStandard(Nsim_lrm::Int=200, mode::FinancialMonteCarlo.SerialMode=FinancialMonteCarlo.SerialMode())
        @assert Nsim_lrm > 0 "Error, number of simulation for VibratoMonteCarloStandard must be positive."
        return new(Nsim_lrm, mode)
    end
end

function init_lrm_vec(x::VibratoMonteCarloStandard, num_typed)
    return Array{typeof(FinancialMonteCarlo.get_rng_type(num_typed))}(undef, x.Nsim_lrm)
end

struct VibratoMonteCarloAnalytic <: AbstractVibrato
    Npoints::Int
    x_min::Float64
    x_max::Float64
    function VibratoMonteCarloAnalytic(Npoints::Int=200, x_min::Float64=-5.0, x_max::Float64=5.0)
        @assert Npoints > 0 "Error, number of integration points for VibratoMonteCarloAnalytic must be positive."
        return new(Npoints, x_min, x_max)
    end
end

function init_lrm_vec(x::VibratoMonteCarloAnalytic, ::Any)
    return range(x.x_min, length=x.Npoints, stop=x.x_max)
end

struct VibratoMonteCarloGaussHermite <: AbstractVibrato
    Npoints::Int
    function VibratoMonteCarloGaussHermite(Npoints::Int=200)
        @assert Npoints > 0 "Error, number of integration points for VibratoMonteCarloGaussHermite must be positive."
        return new(Npoints)
    end
end

function init_lrm_vec(x::VibratoMonteCarloGaussHermite, ::Any)
    return gausshermite(x.Npoints)
end

