module VibratoMonteCarloForwardDiffExt
using VibratoMonteCarlo, ForwardDiff
import VibratoMonteCarlo: v_value

# Recursive, so that nested duals (e.g. ForwardDiff.hessian) lose every derivative level.
v_value(x::ForwardDiff.Dual{Tg, T, N}) where {Tg, T <: Real, N} = v_value(ForwardDiff.value(x))

end # module
