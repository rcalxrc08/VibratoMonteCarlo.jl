module VibratoMonteCarloForwardDiffExt
using VibratoMonteCarlo, ForwardDiff
import VibratoMonteCarlo: v_value

v_value(x::ForwardDiff.Dual{Tg, T, N}) where {Tg, T <: Real, N} = ForwardDiff.value(x)

end # module
