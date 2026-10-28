module VibratoMonteCarloDualNumbersExt
using VibratoMonteCarlo, DualNumbers
import VibratoMonteCarlo: v_value, v_mod

# Recursive for nested duals, which use the generic v_mod.
v_value(x::Dual) = v_value(DualNumbers.value(x))
function v_mod(x::Dual{T}) where {T <: AbstractFloat}
    return dual(one(T), x.epsilon)
end

end # module
