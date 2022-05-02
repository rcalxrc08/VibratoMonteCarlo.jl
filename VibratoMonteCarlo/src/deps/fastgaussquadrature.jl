using .FastGaussQuadrature


function lrm_vec_gauss_analytic!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData, r)
	x,w=Z
	integrand(z) = integrand_lrm(z, mu, sigma, eu_opt, r)
	return sum(@. w*integrand(x*sqrt(2)))/sqrt(pi)
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
