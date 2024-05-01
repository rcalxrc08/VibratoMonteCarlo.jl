using .FastGaussQuadrature

struct VibratoMonteCarloGaussHermite <: AbstractVibrato
    Npoints::Int
    function VibratoMonteCarloGaussHermite(Npoints::Int = 200)
        @assert Npoints > 0 "Error, number of integration points for VibratoMonteCarloGaussHermite must be positive."
        return new(Npoints)
    end
end

function lrm_interface_vec!(mu, sigma, received_args, mcBaseData, mc::VibratoMonteCarloGaussHermite)
    Z = received_args.Z
    eu_opt = received_args.option
    x, w = Z
    args = IntegrandArguments(mu, sigma, eu_opt)
    output = sum(@. integrand_lrm_vec(x * sqrt2, args) * w) / sqrtπ
    # output = sum(w_ * integrand_lrm(x_ * sqrt2, mu, sigma, eu_opt) for (x_, w_) in zip(x, w)) / sqrtπ
    return output

    # # randn!(mcBaseData.parallelMode.rng, Z)
    # #TODO
    # # Z, eu_opt, mcBaseData, mc = received_args
    # Z = received_args.a
    # eu_opt = received_args.b
    # # mcBaseData = received_args.c
    # # mc = received_args.d
    # # @show typeof(sigma)
    # # @show typeof(mu)
    # args = TmpArgs2(mu, sigma, eu_opt)
    # # @show typeof(Z)
    # output = mean(integrand_lrm_vec.(Z, args))
    # return output
end

function init_lrm_vec(x::VibratoMonteCarloGaussHermite, ::Any)
    return gausshermite(x.Npoints)
end
