using .DualNumbers

v_value(x::Dual) = DualNumbers.value(x)
v_mod(x::Dual) = x.epsilon