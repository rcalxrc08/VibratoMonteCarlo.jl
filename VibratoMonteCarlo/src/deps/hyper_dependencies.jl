using .HyperDualNumbers

v_value(x::Hyper) = x.value
v_mod(x::Hyper) = dual(x.epsilon1, x.epsilon12 + x.epsilon1 * x.epsilon2)