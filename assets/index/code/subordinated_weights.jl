# This file was generated, do not modify it. # hide
#hideall
@variables β ψ
mu_sub = log(S_0) + (r - d - ψ) * T + β * δ^2
println("\$\$w_{\\Delta} = ", latex_weight(mu_sub, s_bs, S_0, 1, x, Z), "\$\$")
println("\$\$w_{\\Gamma} = ", latex_weight(mu_sub, s_bs, S_0, 2, x, Z), "\$\$")