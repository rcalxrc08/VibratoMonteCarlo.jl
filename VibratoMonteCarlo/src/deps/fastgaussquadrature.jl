using .FastGaussQuadrature

struct VibratoMonteCarloGaussHermite <: AbstractVibrato
    Npoints::Int
    function VibratoMonteCarloGaussHermite(Npoints::Int=200)
        @assert Npoints > 0 "Error, number of integration points for VibratoMonteCarloGaussHermite must be positive."
        return new(Npoints)
    end
end

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData, r,mc::VibratoMonteCarloGaussHermite)
	x,w=Z
	integrand(z) = integrand_lrm(z, mu, sigma, eu_opt, r)
	return sum(@. w*integrand(x*sqrt(2)))/sqrt(pi)
end

function init_lrm_vec(x::VibratoMonteCarloGaussHermite, ::Any)
    return gausshermite(x.Npoints)
end
