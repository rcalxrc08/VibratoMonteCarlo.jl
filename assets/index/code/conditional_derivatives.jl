# This file was generated, do not modify it. # hide
#hideall
N(x) = 0.0
φ(x) = 0.0
@register_symbolic N(x)
@register_symbolic φ(x)
@register_derivative N(x) 1 φ(x)
@variables S K h [latexwrapper = identity]
@variables D1(S) D2(S)
@variables d_1 [latexwrapper = identity] d_2 [latexwrapper = identity]
# d_1 and d_2 have the same derivative in S.
function generate_delta(price)
    d_prime = 1 / (S * σ * sqrt(h))
    delta = substitute(expand_derivatives(Differential(S)(price)), Dict(Differential(S)(D1) => d_prime, Differential(S)(D2) => d_prime))
    # The identity S e^{-d h} φ(d_1) = K e^{-r h} φ(d_2).
    return simplify(substitute(delta, Dict(φ(D2) => S * exp(-d * h) * φ(D1) / (K * exp(-r * h)))))
end
latex_delta(delta) = latexify(substitute(delta, Dict(D1 => d_1, D2 => d_2)), env = :raw, mult_symbol = "")
delta_call = generate_delta(S * exp(-d * h) * N(D1) - K * exp(-r * h) * N(D2))
delta_binary = substitute(expand_derivatives(Differential(S)(exp(-r * h) * N(D2))), Dict(Differential(S)(D2) => 1 / (S * σ * sqrt(h))))
println("\$\$\\partial_{S} C_{\\text{call}} = ", latex_delta(delta_call), "\$\$")
println("\$\$\\partial_{S} C_{\\text{binary}} = ", latex_delta(delta_binary), "\$\$")