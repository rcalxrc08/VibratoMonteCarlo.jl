using .DualNumbers

v_value(x::Dual) = DualNumbers.value(x)
function v_mod(x::Dual{T}) where {T}
    return dual(one(T), x.epsilon)
end