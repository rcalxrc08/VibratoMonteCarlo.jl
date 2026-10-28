module VibratoMonteCarloFastGaussQuadratureExt
using VibratoMonteCarlo, FastGaussQuadrature
import VibratoMonteCarlo: init_lrm_vec, VibratoMonteCarloGaussHermite

function init_lrm_vec(x::VibratoMonteCarloGaussHermite, ::Any)
    return gausshermite(x.Npoints)
end

end # module
