# This file was generated, do not modify it. # hide
#hideall
@variables λ κ
with_jumps(delta) = substitute(delta, Dict(d => d + λ * κ))
println("\$\$\\partial_{S} C_{\\lambda \\kappa,\\text{call}} = ", latex_delta(with_jumps(delta_call)), "\$\$")
println("\$\$\\partial_{S} C_{\\lambda \\kappa,\\text{binary}} = ", latex_delta(with_jumps(delta_binary)), "\$\$")