module VibratoMonteCarloMooncakeExt
using VibratoMonteCarlo, Mooncake
import VibratoMonteCarlo: v_value, v_mod
using Mooncake: @mooncake_overlay, @zero_derivative, MinimalCtx

# Mooncake differentiates the floats themselves: v_value has to stop the gradient, and v_mod needs the exp form,
# since v_mod(::AbstractFloat) = 1 has no derivative and drops every score term.
@zero_derivative MinimalCtx Tuple{typeof(v_value), Base.IEEEFloat}
@mooncake_overlay v_mod(x::Base.IEEEFloat) = exp(x - v_value(x))

# The pivoted QR of the regressions of the exercise policy ends in a foreigncall, and it only sees values.
@zero_derivative MinimalCtx Tuple{typeof(VibratoMonteCarlo.exercise_regions), Vararg}

end # module
