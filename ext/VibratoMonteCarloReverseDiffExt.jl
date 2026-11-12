module VibratoMonteCarloReverseDiffExt
using VibratoMonteCarlo, ReverseDiff
import VibratoMonteCarlo: v_value

function v_value(x::ReverseDiff.TrackedReal{V, D, O}) where {V <: Real, D <: Real, O}
    return v_value(x.value)
end

end # module
