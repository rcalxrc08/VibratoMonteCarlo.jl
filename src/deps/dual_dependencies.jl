using .DualNumbers

v_value(x::Dual) = DualNumbers.value(x)
v_mod(x::Dual) = dual(1.0, x.epsilon)