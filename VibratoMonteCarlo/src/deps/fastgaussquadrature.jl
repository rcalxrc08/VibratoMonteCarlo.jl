using .FastGaussQuadrature

struct VibratoMonteCarloGaussHermite <: AbstractVibrato
    Npoints::Int
    function VibratoMonteCarloGaussHermite(Npoints::Int = 200)
        @assert Npoints > 0 "Error, number of integration points for VibratoMonteCarloGaussHermite must be positive."
        return new(Npoints)
    end
end

function lrm_interface!(Z, mu, sigma, eu_opt::FinancialMonteCarlo.EuropeanPayoff, mcBaseData, r, mc::VibratoMonteCarloGaussHermite)
    x, w = Z
    zero_typed = FinancialMonteCarlo.predict_output_type_zero(mu, sigma, eu_opt)
    integrand(z)::typeof(zero_typed) = integrand_lrm(z, mu, sigma, eu_opt, r)
    output = sum(w_ * integrand(x_ * sqrt(2)) for (x_, w_) in zip(x, w)) / sqrt(pi)
    return output
end

function init_lrm_vec(x::VibratoMonteCarloGaussHermite, ::Any)
    return gausshermite(x.Npoints)
end
