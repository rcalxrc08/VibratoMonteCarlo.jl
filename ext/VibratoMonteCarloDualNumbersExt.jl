module VibratoMonteCarloDualNumbersExt
using VibratoMonteCarlo, DualNumbers
import VibratoMonteCarlo: v_value, v_mod

v_value(x::Dual) = DualNumbers.value(x)
function v_mod(x::Dual{T}) where {T}
    return dual(one(T), x.epsilon)
end

end # module
