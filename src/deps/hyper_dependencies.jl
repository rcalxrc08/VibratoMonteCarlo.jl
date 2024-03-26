using .HyperDualNumbers

v_value(x::Hyper) = x.value
function v_mod(x::Hyper{T}) where {T}
    der_2 = @muladd x.epsilon12 + x.epsilon1 * x.epsilon2
    return hyper(one(T), x.epsilon1, x.epsilon2, der_2)
end