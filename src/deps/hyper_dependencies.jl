using .HyperDualNumbers

v_value(x::Hyper) = x.value
v_mod(x::Hyper) = hyper(1.0, x.epsilon1, x.epsilon2, x.epsilon12 + x.epsilon1 * x.epsilon2)