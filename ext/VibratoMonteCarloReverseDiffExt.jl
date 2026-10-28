module VibratoMonteCarloReverseDiffExt
using VibratoMonteCarlo, ReverseDiff
import VibratoMonteCarlo: v_value

function v_value(x::ReverseDiff.TrackedReal{V, D, O}) where {V <: Real, D <: Real, O}
    return v_value(x.value)
end

# using ChainRulesCore;
# function ChainRulesCore.rrule(::typeof(VibratoMonteCarlo.v_value), x::AbstractFloat)
    # y = VibratoMonteCarlo.v_value(x)
    # function v_value_pullback(ȳ)
        # return ChainRulesCore.NoTangent(), zero(ȳ)
    # end
    # return y, v_value_pullback
# end

# ReverseDiff.@grad_from_chainrules VibratoMonteCarlo.v_value(x)

# function ChainRulesCore.rrule(::typeof(VibratoMonteCarlo.v_mod), x::AbstractFloat)
    # y = VibratoMonteCarlo.v_mod(x)
    # function v_mod_pullback(ȳ)
        # return ChainRulesCore.NoTangent(), ȳ
    # end
    # return y, v_mod_pullback
# end

# ReverseDiff.@grad_from_chainrules VibratoMonteCarlo.v_mod(x)

end # module
